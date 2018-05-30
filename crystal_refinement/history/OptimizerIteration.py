from __future__ import absolute_import
import copy, random
from graphviz import Digraph
from collections import OrderedDict
from itertools import chain
from crystal_refinement.utils import optimizer_scores as scores


class OptimizerIteration:

    """
    Define class to hold information on one optimizer iteration
    """
    def __init__(self, parent, ins_file, res_file, r1, bond_score, n_missing_elements, stoich_score, anisotropy_penalty,
                 overall_score, score_weighting=1.0, annotation=None):
        self.ins_file = copy.deepcopy(ins_file)
        self.res_file = copy.deepcopy(res_file)

        self.parent = parent
        self.dead_branch = False
        self.children = []
        self.annotation = annotation

        # the higher this is, the more the r1 counts. Must be between 0 and 1
        self.score_weighting = score_weighting
        self.r1 = r1
        self.bond_score = bond_score
        self.n_missing_elements = n_missing_elements
        self.stoich_score = stoich_score
        self.anisotropy_penalty = anisotropy_penalty
        self.overall_score = overall_score

    @classmethod
    def build_with_bond_list(cls, parent, ins_file, res_file, bonds, cache, score_weighting=1.0, annotation=None):
        r1 = res_file.r1
        bond_score = scores.get_compound_bond_score(bonds, res_file, cache)
        n_missing_elements = scores.get_missing_element_score(res_file, cache)
        stoich_score = scores.get_stoichiometry_score(res_file, cache)
        anisotropy_penalty = scores.get_negative_anisotropic_penalty(res_file)
        return cls(
            parent=parent,
            ins_file=copy.deepcopy(ins_file),
            res_file=copy.deepcopy(res_file),
            r1=r1,
            bond_score=bond_score,
            n_missing_elements=n_missing_elements,
            stoich_score=stoich_score,
            anisotropy_penalty=anisotropy_penalty,
            overall_score=scores.get_overall_score(r1,
                                                   bond_score,
                                                   n_missing_elements,
                                                   stoich_score,
                                                   anisotropy_penalty,
                                                   score_weighting),
            score_weighting=score_weighting,
            annotation=annotation
        )

    def add_child(self, child):
        """
        Add a child iteration to the current iteration node
        :param child: iteration to add as a child
        """
        self.children.append(child)

    def get_full_history(self):
        """
        Get the path of nodes from the root node to this node
        """
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
        """
        Get all leaf nodes with this node as the root.
        :return:
        """
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
        """
        Get a copy of the res file associated with this node.
        :return:
        """
        return copy.deepcopy(self.res_file)

    def get_ins_copy(self):
        """
        Get a copy of the ins file associated with this node.
        :return:
        """
        return copy.deepcopy(self.ins_file)

    def generate_graph(self, output_file):
        """
        Generate a graphviz image of the tree with this node as the root
        :param output_file: the file destination to output the image
        :return:
        """
        dot = Digraph()
        node_label = str(random.getrandbits(32))
        dot.node(node_label, str(self.r1), color="green")
        sorted_leaves = self.get_sorted_leaves()
        for i, child in enumerate(self.children):
            child._generate_graph(str(random.getrandbits(32)), node_label, dot, sorted_leaves)
        dot.render(output_file, view=False)

    def _generate_graph(self, node_label, parent_label, dot, sorted_leaves):
        """
        Actually do the work to generate the graph. This function is recursive.
        :param node_label: graphviz label for the current node. This is supposed to be unique among the nodes in the
        graph, so we just generate a random one.
        :param parent_label: graphviz label for the parent node so that we can hook them up
        :param dot: the graphviz object
        :param sorted_leaves: ranking of the leaves so we can label nodes with rank when necessary.
        :return:
        """
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
            dot.node(node_label, self.generate_label(rank=rank+1), color=color)
        else:
            dot.node(node_label, self.generate_label(), color=color)
        dot.edge(parent_label, node_label, color=color, label=self.annotation)

        return highlight

    def generate_truncated_graph(self, output_file):
        """
        Generate a graphviz image of a truncated tree with this node as the root. If all children of a node are dead
        branches, the children for that node will not be rendered.
        :param output_file: the file destination to output the image
        :return:
        """
        dot = Digraph()
        node_label = str(random.getrandbits(32))
        dot.node(node_label, self.generate_label(), color="green")
        sorted_leaves = self.get_sorted_leaves()
        self._generate_truncated_graph(node_label, dot, sorted_leaves)
        dot.render(output_file, view=False)

    def _generate_truncated_graph(self, node_label, dot, sorted_leaves):
        """
        Actually do the work to generate the truncated graph. This function is recursive.
        :param node_label: graphviz label for the current node. This is supposed to be unique among the nodes in the
        graph, so we just generate a random one.
        :param dot: the graphviz object
        :param sorted_leaves: ranking of the leaves so we can label nodes with rank when necessary.
        :return:
        """

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
                    dot.node(child_label, self.generate_label(rank=rank + 1), color=color)
                except ValueError:
                    dot.node(child_label, self.generate_label(), color=color)
            else:
                dot.node(child_label, self.generate_label(), color=color)
            dot.edge(node_label, child_label, color=color, label=child.annotation)

        if len(self.children) == 0:
            rank = sorted_leaves.index(self)
            if rank == 0:
                highlight = True

        return highlight

    def is_dead_branch(self):
        """
        :return: If this node is in a dead branch
        """
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
        """
        Update the dead branch labels of all nodes under this root. If a node has no children, it is labeled as dead.
        :return:
        """
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

    def propagate(self):
        """
        Copy this iteration into the next generation.
        :return:
        """
        new_annotation = None
        if self.annotation is not None:
            new_annotation = "Propagated from previous generation"
        new_child = OptimizerIteration(self,
                                       self.ins_file,
                                       self.res_file,
                                       self.r1,
                                       self.bond_score,
                                       self.n_missing_elements,
                                       self.stoich_score,
                                       self.anisotropy_penalty,
                                       self.overall_score,
                                       self.score_weighting,
                                       new_annotation)

        self.children.append(new_child)

    def get_score(self, criterion="overall_score"):
        """
        Get the score associated with the .res file of this node
        :param criterion:
        :return:
        """
        assert criterion in ["overall_score", "r1_only", "bond_only", "missing_elements", "stoich_only", "anisotropy_only"],\
            "{} is not a supported criterion".format(criterion)
        if criterion == "overall_score":
            return self.overall_score
        if criterion == "r1_only":
            return self.r1
        if criterion == "bond_only":
            return self.bond_score
        if criterion == "missing_elements":
            return self.n_missing_elements
        if criterion == "stoich_only":
            return self.stoich_score
        if criterion == "anisotropy_only":
            return self.anisotropy_penalty

    def get_sorted_leaves(self, criteria="overall_score"):
        """
        Return the leaf nodes underneath this node sorted by the provided criteria.
        :param criteria: to sort by
        :return:
        """
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

    def generate_label(self, rank=None):
        """
        Generate the label for this node in the graph.
        :param rank: Optional rank
        :return: label
        """
        rank_label = ""
        if rank is not None:
            rank_label = "\nrank:{}".format(rank)
        penalty_report = ""
        if self.n_missing_elements > 0:
            penalty_report += "\nmissing_elements: {}".format(self.n_missing_elements)
        if self.anisotropy_penalty > 0:
            penalty_report += "\nnegative anisotropy: {}".format(self.anisotropy_penalty)
        if len(penalty_report) > 0:
            penalty_report = "\npenalties" + penalty_report
        return "r1: {}\nbond: {}\nstochiometry: {}{}\noverall:{}{}".format(
            self.r1, self.bond_score, self.stoich_score, penalty_report, self.get_score(), rank_label)

    def get_best(self, criteria="overall_score"):
        """
        Return the best leaf node under this one according to the provided critera
        :return:
        """
        return self.get_sorted_leaves(criteria)[0]