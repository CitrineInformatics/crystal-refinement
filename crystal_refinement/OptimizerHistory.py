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
        self.overall_score = self.r1 ** 3 * self.bond_score * 1e4
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
        node_label = str(random.getrandbits(32))
        dot.node(node_label, "r1={}, bonds={:.4}\noverall={:.5}".format(self.r1, self.bond_score, self.overall_score), color="green")
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
        if highlight:
            dot.node(node_label, "r1={}, bonds={:.4}\noverall={:.5}".format(self.r1, self.bond_score, self.overall_score), color="green")
            dot.edge(parent_label, node_label, color="green", label=self.annotation)
        else:
            dot.node(node_label, "r1={} bonds={:.4}\noverall={:.5}".format(self.r1, self.bond_score, self.overall_score))
            dot.edge(parent_label, node_label, label=self.annotation)
        return highlight


    # Copy this iteration into the next generation
    def propagate(self):
        new_annotation = None
        if self.annotation is not None:
            new_annotation = "Propagated from previous generation"
        new_child = OptimizerIteration(self, self.get_ins(), self.get_res(), self.bond_score, new_annotation)
        self.children.append(new_child)

    def get_best(self):
        return sorted(self.get_leaves(), key=lambda iteration: iteration.overall_score)[0]


class OptimizerHistory:
    """
    Define class to hold information optimizer history information
    """
    def __init__(self, driver, ins_file):
        self.driver = driver
        res = self.driver.run_SHELXTL(ins_file)
        bonds = utils.get_bonds(self.driver, res)
        self.head = OptimizerIteration(None, ins_file, res, utils.score_compound_bonds(bonds, ins_file))
        self.leaves = [self.head]

    def run_iter(self, ins_file, parent_iteration, annotation=None):
        """
        Run the given ins file through SHELXTL and record the file and resulting r1

        :param ins_file: SHELXFile object
        :return res: SHELXFile object
        """
        res = self.driver.run_SHELXTL(ins_file)
        if res is None:
            return None
        bonds = utils.get_bonds(self.driver, res)
        new_iter = OptimizerIteration(parent_iteration, ins_file, res, utils.score_compound_bonds(bonds, ins_file), annotation)
        return new_iter

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
