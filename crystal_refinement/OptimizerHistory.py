import copy, random, utils
from graphviz import Digraph
import numpy as np


class OptimizerIteration:
    """
    Define class to hold information on one optimizer iteration
    """
    def __init__(self, parent, ins_file, res_file, bond_score, annotation=None):
        self.ins_file = copy.deepcopy(ins_file)
        self.res_file = copy.deepcopy(res_file)
        self.r1 = res_file.r1
        self.bond_score = bond_score
        self.overall_score = self.r1 * self.bond_score
        self.parent = parent
        self.dead_branch = False
        self.children = []
        self.annotation = annotation

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

    def get_res(self):
        return copy.deepcopy(self.res_file)

    def get_ins(self):
        return copy.deepcopy(self.ins_file)

    def generate_graph(self, output_file):
        dot = Digraph()
        node_label = str(random.getrandbits(16))
        dot.node(node_label, str(self.r1), color="green")
        best = self.get_best()
        for i, child in enumerate(self.children):
            child._generate_graph(str(random.getrandbits(16)), node_label, dot, best)
        dot.render(output_file, view=False)

    def _generate_graph(self, node_label, parent_label, dot, best):
        highlight = False
        for i, child in enumerate(self.children):
            if child._generate_graph(str(random.getrandbits(16)), node_label, dot, best):
                highlight = True
        if len(self.children) == 0:
            if self == best:
                highlight = True
        color = "black"
        if highlight:
            color = "green"
        if self.dead_branch:
            color = "red"

        dot.node(node_label, str(self.r1), color=color)
        dot.edge(parent_label, node_label, color=color, label=self.annotation)

        return highlight


    # Copy this iteration into the next generation
    def propagate(self):
        new_annotation = None
        if self.annotation is not None:
            new_annotation = "Propagated from previous generation"
        new_child = OptimizerIteration(self, self.get_ins(), self.get_res(), self.bond_score, new_annotation)
        self.children.append(new_child)

    def get_sorted_leaves(self):
        return sorted(self.get_leaves(), key=lambda iteration: (iteration.r1, len(iteration.res_file.mixed_site_numbers)))

    def get_best(self):
        return self.get_sorted_leaves()[0]


class OptimizerHistory:
    """
    Define class to hold information optimizer history information
    """
    def __init__(self, driver, ins_file, max_n_leaves=10):
        self.driver = driver
        res = self.driver.run_SHELXTL(ins_file)
        bonds = utils.get_bonds(self.driver, res)
        self.head = OptimizerIteration(None, ins_file, res, utils.score_compound_bonds(bonds, ins_file))
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
            bonds = utils.get_bonds(self.driver, res)
        except IndexError:
            return None
        new_iter = OptimizerIteration(parent_iteration, ins_file, res, utils.score_compound_bonds(bonds, ins_file), annotation)
        return new_iter

    def clean_history(self):
        if len(self.leaves) > self.max_n_leaves:
            sorted_leaves = self.head.get_sorted_leaves()
            for i in range(self.max_n_leaves, len(self.leaves)):
                # print "Dropping leaf: "
                # print sorted_leaves[i].ins_file.get_crystal_sites_text()
                sorted_leaves[i].dead_branch = True
            self.update_leaves()


    def save(self, iterations):
        if not isinstance(iterations, list):
            iterations = [iterations]
        for iteration in iterations:
            iteration.parent.add_child(iteration)
        self.update_leaves()

    def run_and_save(self, ins_file, parent_iteration):
        iteration = self.run_iter(ins_file, parent_iteration)
        if iteration is not None:
            self.save(iteration)
        return iteration

    def update_leaves(self):
        self.leaves = self.head.get_leaves()

    def get_best_history(self):
        return self.head.get_best().get_full_history()
