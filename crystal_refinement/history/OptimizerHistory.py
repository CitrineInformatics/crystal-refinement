from crystal_refinement.history.OptimizerIteration import OptimizerIteration
from crystal_refinement.utils.bond_utils import get_bonds


class OptimizerHistory:
    """
    Define class to hold information optimizer history information
    """
    def __init__(self, driver, cache, ins_file, score_weighting=1.0, max_n_leaves=50):
        self.driver = driver
        self.cache = cache
        self.score_weighting = score_weighting
        res = self.driver.run_SHELXTL(ins_file)
        bonds = get_bonds(self.driver, res)
        self.head = OptimizerIteration.build_with_bond_list(None, ins_file, res, bonds, self.cache,
                                                            score_weighting=self.score_weighting)

        self.max_n_leaves = max_n_leaves

    def run_iter(self, ins_file, parent_iteration, annotation=None):
        """
        Run the given ins file through SHELXTL and record the file and resulting r1

        :param ins_file: SHELXFile object
        :param parent_iteration: the parent OptimizerIteration object
        :param annotation: the annotation for the new iteration
        :return res: SHELXFile object
        """
        res = self.driver.run_SHELXTL(ins_file)
        if res is None:
            return None
        # If refinement is unstable, no cif file will be generated. The iteration should fail then
        try:
            bonds = get_bonds(self.driver, res)
            if len(bonds) == 0:
                return None
            new_iter = OptimizerIteration.build_with_bond_list(parent_iteration,
                                                              ins_file,
                                                              res,
                                                              bonds,
                                                              self.cache,
                                                              score_weighting=self.score_weighting,
                                                              annotation=annotation)
        # This throws a FileNotFound error in python3, let's just try catching everything for now...
        except:
            return None

        return new_iter

    def clean_history(self, n_to_keep=None, branch=None, criteria="overall_score"):
        """
        Prune the tree according to overall scores.
        :param n_to_keep: Number of leaves to keep.
        :param branch: Branch to prune.
        :param criteria: Criteria with which to prune the branches.
        """
        if branch is None:
            branch = self.head
        if n_to_keep is None:
            n_to_keep = self.max_n_leaves
        if len(branch.get_leaves()) > n_to_keep:
            sorted_leaves = branch.get_sorted_leaves(criteria)
            for i in range(n_to_keep, len(branch.get_leaves())):
                sorted_leaves[i].dead_branch = True

    def save(self, iterations):
        """
        Save the specified iterations in the optimizer history tree
        :param iterations: iterations to save
        :return:
        """
        if not isinstance(iterations, list):
            iterations = [iterations]
        for iteration in iterations:
            iteration.parent.add_child(iteration)

    def run_and_save(self, ins_file, parent_iteration, annotation=None):
        """
        Run SHELX on the specified ins file and save the result under the specified parent iteration. If there is an
        issue with running SHELX (iteration == None), don't save.
        :param ins_file:
        :param parent_iteration:
        :param annotation:
        :return:
        """
        iteration = self.run_iter(ins_file, parent_iteration, annotation)
        if iteration is not None:
            self.save(iteration)
        return iteration

    def get_leaves(self):
        """
        :return: the leaves for the full history tree
        """
        return self.head.get_leaves()

    def get_best_history(self):
        """
        :return: The history of the best leaf (according to the overall score criteria)
        """
        return self.head.get_best().get_full_history()

    def generate_graph(self, output_file):
        """
        Generate the graph representation of the history tree.
        :param output_file: file where the graph will be output
        :return:
        """
        self.head.generate_graph(output_file)
