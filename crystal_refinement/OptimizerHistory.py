import copy, random, math
from graphviz import Digraph
from collections import OrderedDict
from itertools import chain


class OptimizerIteration:
    """
    Define class to hold information on one optimizer iteration
    """
    def __init__(self, parent, ins_file, res_file, bond_score, n_missing_elements, stoich_score, score_weighting=1.0, annotation=None):
        self.ins_file = copy.deepcopy(ins_file)
        self.res_file = copy.deepcopy(res_file)
        self.r1 = res_file.r1
        self.bond_score = bond_score
        # the higher this is, the more the r1 counts. Must be between 0 and 1
        self.score_weighting = score_weighting

        self.parent = parent
        self.dead_branch = False
        self.children = []
        self.annotation = annotation
        self.n_missing_elements = n_missing_elements
        self.stoich_score = stoich_score

    def add_child(self, child):
        self.children.append(child)

    def get_full_history(self):
        history = []
        cur_iter = self
        while True:
            history.append(cur_iter)
            cur_iter = cur_iter.parent
            if cur_iter is None:
                break
        history.reverse()
        return history

    def get_leaves(self):
        leaves = []
        if len(self.children) == 0:
            if self.dead_branch:
                return []
            else:
                return [self]

        for child in self.children:
            leaves.extend(child.get_leaves())
        return leaves

    def get_res_copy(self):
        return copy.deepcopy(self.res_file)

    def get_ins_copy(self):
        return copy.deepcopy(self.ins_file)

    def generate_graph(self, output_file):
        dot = Digraph()
        node_label = str(random.getrandbits(32))
        dot.node(node_label, str(self.r1), color="green")
        sorted_leaves = self.get_sorted_leaves()
        for i, child in enumerate(self.children):
            child._generate_graph(str(random.getrandbits(32)), node_label, dot, sorted_leaves)
        dot.render(output_file, view=False)

    def _generate_graph(self, node_label, parent_label, dot, sorted_leaves):
        highlight = False
        for i, child in enumerate(self.children):
            if child._generate_graph(str(random.getrandbits(32)), node_label, dot, sorted_leaves):
                highlight = True
        if len(self.children) == 0:
            try:
                rank = sorted_leaves.index(self)
            except ValueError:
                rank = -1
            if rank == 0:
                highlight = True
        color = "black"
        if highlight:
            color = "green"
        if self.dead_branch:
            color = "red"
        if len(self.children) == 0 and rank >= 0:
            dot.node(node_label, "r1: {}\nbond: {}\nmissing_elements: {}\noverall:{}\nrank:{}"
                     .format(self.r1, self.bond_score, self.n_missing_elements, self.get_score(), rank + 1), color=color)
        else:
            dot.node(node_label, "r1: {}\nbond: {}\nmissing_elements: {}\noverall:{}".format(self.r1, self.bond_score, self.n_missing_elements, self.get_score()),
                     color=color)
        dot.edge(parent_label, node_label, color=color, label=self.annotation)

        return highlight


    def generate_truncated_graph(self, output_file):
        dot = Digraph()
        node_label = str(random.getrandbits(32))
        dot.node(node_label, "r1: {}\nbond: {}\nmissing_elements: {}\noverall:{}".format(self.r1, self.bond_score, self.n_missing_elements, self.get_score()), color="green")
        sorted_leaves = self.get_sorted_leaves()
        self._generate_truncated_graph(node_label, dot, sorted_leaves)
        dot.render(output_file, view=False)

    def _generate_truncated_graph(self, node_label, dot, sorted_leaves):
        # The parent node (this) controls the rendering of it's child nodes
        # if all children are dead branches, don't render any children
        highlight = False
        for child in self.children:
            child_label = str(random.getrandbits(32))
            color = "black"
            if child.is_dead_branch():
                color = "red"
            else:
                if child._generate_truncated_graph(child_label, dot, sorted_leaves):
                    color = "green"
                    highlight = True
            if len(child.children) == 0:
                try:
                    rank = sorted_leaves.index(child)
                    dot.node(child_label, "r1: {}\nbond: {}\nmissing_elements: {}\noverall:{}\nrank:{}"
                     .format(self.r1, self.bond_score, self.n_missing_elements, self.get_score(), rank + 1),
                             color=color)
                except ValueError:
                    dot.node(child_label,
                             "r1: {}\nbond: {}\nmissing_elements: {}\noverall:{}".format(self.r1, self.bond_score, self.n_missing_elements, self.get_score()),
                             color=color)
            else:
                dot.node(child_label, "r1: {}\nbond: {}\nmissing_elements: {}\noverall:{}".format(self.r1, self.bond_score, self.n_missing_elements, self.get_score()),
                         color=color)
            dot.edge(node_label, child_label, color=color, label=child.annotation)

        if len(self.children) == 0:
            rank = sorted_leaves.index(self)
            if rank == 0:
                highlight = True

        return highlight

    def is_dead_branch(self):
        if self.dead_branch:
            return True
        else:
            if len(self.children) == 0:
                return False
            for child in self.children:
                if not child.is_dead_branch():
                    return False
        return True

    def update_dead_branches(self):

        if self.is_dead_branch():
            self.dead_branch = True
        else:
            if len(self.children) == 0:
                self.dead_branch = False
            else:
                dead_branch = True
                for child in self.children:
                    child.update_dead_branches()
                    if not child.is_dead_branch():
                        dead_branch = False
                self.dead_branch = dead_branch

    # Copy this iteration into the next generation
    def propagate(self):
        new_annotation = None
        if self.annotation is not None:
            new_annotation = "Propagated from previous generation"
        new_child = OptimizerIteration(self, self.get_ins_copy(), self.get_res_copy(), self.bond_score, self.n_missing_elements, self.stoich_score, self.score_weighting,
        new_annotation)
        self.children.append(new_child)

    def get_score(self, criterion="overall_score"):
        assert criterion in ["overall_score", "r1_only", "bond_only", "missing_elements"], "{} is not a supported criterion".format(
            criterion)
        if criterion == "overall_score":
            bond_score_basis = 0.1
            missing_elements_basis = 0.05
            stoich_score_basis = 0.1
            # return math.pow(self.r1, self.score_weighting) + math.pow(self.bond_score, 1 - self.score_weighting)
            # return (1.0 + math.pow(self.r1 / r1_basis, self.score_weighting)) \
            #        * (1.0 + math.pow(self.bond_score / bond_score_basis, 1 - self.score_weighting))
            # try:
            return math.pow(self.r1, self.score_weighting) * \
                   math.pow(self.bond_score + bond_score_basis, 1 - self.score_weighting) + \
                   missing_elements_basis * self.n_missing_elements ** 2 + \
                   stoich_score_basis * self.stoich_score
            # except TypeError:
            #     print(self.r1, self.score_weighting, self.bond_score, self.n_missing_elements)
            #     quit()
            # return self.r1 * self.bond_score
            # return (self.r1, len(self.res_file.mixed_site_numbers))
        if criterion == "r1_only":
            return self.r1
        if criterion == "bond_only":
            return self.bond_score
        if criterion == "missing_elements":
            return self.n_missing_elements
        if criterion == "stoich_only":
            return self.stoich_score

    def get_sorted_leaves(self, criteria="overall_score"):
        if type(criteria) != list:
            criteria = [criteria]

        sorted_by_criteria = []

        for criterion in criteria:
            sorted_by_criteria.append(sorted(self.get_leaves(), key=lambda iteration: iteration.get_score(criterion)))

        # interleave sorted lists
        sorted_leaves_with_duplicates = list(chain(*zip(*sorted_by_criteria)))

        # deduplicate interleaved list
        deduplicated = list(OrderedDict.fromkeys(sorted_leaves_with_duplicates))

        return deduplicated

    def get_best(self):
        return self.get_sorted_leaves()[0]


class OptimizerHistory:
    """
    Define class to hold information optimizer history information
    """
    def __init__(self, driver, utils, ins_file, score_weighting=1.0, max_n_leaves=50):
        self.driver = driver
        self.utils = utils
        self.score_weighting = score_weighting
        res = self.driver.run_SHELXTL(ins_file)
        bonds = self.utils.get_bonds(self.driver, res)
        self.head = OptimizerIteration(None, ins_file, res, self.utils.score_compound_bonds_relative_distance(bonds, ins_file, self.driver),
                                       self.utils.get_n_missing_elements(res), self.utils.get_stoichiometry_score(res), score_weighting=self.score_weighting)

        self.leaves = [self.head]
        self.max_n_leaves = max_n_leaves

    def run_iter(self, ins_file, parent_iteration, annotation=None):
        """
        Run the given ins file through SHELXTL and record the file and resulting r1

        :param ins_file: SHELXFile object
        :return res: SHELXFile object
        """
        res = self.driver.run_SHELXTL(ins_file)
        if res is None:
            return None
        # If refinement is unstable, no cif file will be generated. The iteration should fail then
        try:
            bonds = self.utils.get_bonds(self.driver, res)
            if len(bonds) == 0:
                return None

            new_iter = OptimizerIteration(parent_iteration, ins_file, res,
                                          self.utils.score_compound_bonds_relative_distance(bonds, res, self.driver),
                                          self.utils.get_n_missing_elements(res),
                                          self.utils.get_stoichiometry_score(res), score_weighting=self.score_weighting,
                                          annotation=annotation)
        # This throws a FileNotFound error in python3, let's just try catching everything for now...
        # except (IndexError, ZeroDivisionError):
        except:
            return None

        return new_iter

    def clean_history(self, n_to_keep=None, branch=None, criteria="overall_score"):
        """
        Prune the tree according to overall scores.
        :param n_to_keep: Number of leaves to keep.
        :param branch: Branch to prune.
        """
        if branch is None:
            branch = self.head
        if n_to_keep is None:
            n_to_keep = self.max_n_leaves
        if len(branch.get_leaves()) > n_to_keep:
            sorted_leaves = branch.get_sorted_leaves(criteria)
            for i in range(n_to_keep, len(branch.get_leaves())):
                # print("Dropping leaf: ")
                # print(sorted_leaves[i].ins_file.get_crystal_sites_text())
                sorted_leaves[i].dead_branch = True
            self.update_leaves()


    def save(self, iterations):
        if not isinstance(iterations, list):
            iterations = [iterations]
        for iteration in iterations:
            iteration.parent.add_child(iteration)
        self.update_leaves()

    def run_and_save(self, ins_file, parent_iteration, annotation=None):
        iteration = self.run_iter(ins_file, parent_iteration, annotation)
        if iteration is not None:
            self.save(iteration)
        return iteration

    def update_leaves(self):
        self.leaves = self.head.get_leaves()

    def get_best_history(self):
        return self.head.get_best().get_full_history()

    def generate_graph(self, output_file):
        self.head.generate_graph(output_file)
