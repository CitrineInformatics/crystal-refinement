from SHELXDriver import SHELXDriver
import os, re, time
import shutil
from OptimizerHistory import OptimizerHistory
from OptimizerSteps import OptimizerSteps
from OptimizerUtils import OptimizerUtils
from collections import defaultdict


class Optimizer:
    """
    Class for performing single crystal refinement
    """
    def __init__(self, path_to_xl, path_to_xs, path_to_ins, input_prefix, output_prefix, use_wine=False,
                 bond_lengths=None, mixing_pairs=None, use_ml_model=False,
                 r1_similarity_threshold=0.0075, occupancy_threshold=0.02, r1_threshold=0.1, score_weighting=0.8,
                 max_n_leaves=50, least_squares_iterations=4, n_results=10):
        """
        :param path_to_xl: path to xl executable (including executable file name)
        :param path_to_xs: path to xs executable (including executable file name)
        :param path_to_ins: path to directory containing *.ins file output by xprep (not including *.ins file name)
        :param input_prefix: prefix of ins file (e.g. for file.ins, the prefix would be "file")
        :param output_prefix: prefix of the result file that the optimizer will output
        :param use_wine: In order to run Windows executables on a mac, the wine app can be used.
        :param bond_lengths: A list of prescribed ideal bond lengths in Angstroms.  Should be in the form of triples.
            (e.g. [('Ag', 'Ag', 2.85)] ).  If a bond length is not specified, then it is calculated via either a machine
            learning model or based on the atomic radii.
        :param mixing_pairs: A list of acceptable elements to allow to co-occupy a site. (e.g. [('Ge', 'Ru')])
            If no such list is provided, then acceptable mixing pairs are determined based on pymatgen information.
        :param use_ml_model: Whether to use the bond length machine learning model to estimate bond lengths.
            If True, a machine learning model hosted at citrination.com is queried to determine the bond lengths.
                Using this option requires obtaining a free Citrination user account and setting the CITRINATION_API_KEY
                environment variable with your account API key.
                If this machine learning model has high uncertainty for the bond lengths queried, the optimizer defaults
                back to the atomic radii model.
            If False, a simpler bond length model based on atomic radii is used.
        :param r1_similarity_threshold: If r1 scores for multiple options are within this similarity threshold, the
            optimizer will branch and explore all the options.  (Should be doulbe in range (0.0, 1.0)).
        :param occupancy_threshold: Minimum deviation in occupancy or site mixing to trigger partial occupancy
            (Should be double in range (0.0, 1.0)).
        :param r1_threshold: Threshold to denote a high r1 score, which triggers an alternate path in the optimizer
            where bond lengths drive the optimization instead of r1 (Should be a double in range (0.0, 1.0)).
        :param score_weighting: Weighting of R1 score versus bond length score when choosing optimal path.
            A value of 1.0 corresponds to R1 only.  (Should be a double in range (0.0, 1.0)).
        :param max_n_leaves: Maximum number of optimization paths to follow at any given optimization step.
            A higher value may give a more optimal path, but will be more computationally expensive.
            (Should be positive integer).
        :param least_squares_iterations: Number of least squares iterations for xl to use
            (Should be a positive integer).
        :param n_results: Number of final .res files to save.  For n_results=1, only the results file with the
            best final score (based on R1 and bond lengths) will be saved.  (Should be a non-negative integer).

        """

        # Initialize parameters based on arguments
        self.path_to_xl = path_to_xl
        self.path_to_xs = path_to_xs
        self.path_to_ins = path_to_ins
        self.input_prefix = input_prefix
        self.output_prefix = output_prefix
        self.use_wine = use_wine
        self.bond_lengths = bond_lengths
        self.mixing_pairs = mixing_pairs
        self.use_ml_model = use_ml_model
        self.overall_score_similarity_threshold = 0.05
        self.r1_similarity_threshold = r1_similarity_threshold
        self.r1_threshold = r1_threshold
        self.occupancy_threshold = occupancy_threshold
        self.least_squares_iterations = least_squares_iterations
        self.score_weighting = score_weighting
        self.max_n_leaves = max_n_leaves
        self.optimizer_steps = OptimizerSteps(self)
        self.n_results = n_results

        self.check_inputs()

        # Initialize objects to be used during run()
        self.driver = SHELXDriver(ins_path=self.path_to_ins, prefix=self.output_prefix, path_to_xl=self.path_to_xl, path_to_xs=self.path_to_xs, use_wine=self.use_wine)
        self.history = None
        self.utils = None

    def run(self):
        """
        Method to run the optimization

        :return:
        """

        # Copy ins and hkl file to output prefix
        os.chdir(self.path_to_ins)
        shutil.copy(os.path.join(self.path_to_ins, self.input_prefix + ".hkl"), os.path.join(self.path_to_ins, self.output_prefix + ".hkl"))
        shutil.copy(os.path.join(self.path_to_ins, self.input_prefix + ".ins"), os.path.join(self.path_to_ins, self.output_prefix + ".ins"))

        # Check that the ins file is direct from xprep, without having been run before
        f = open(self.output_prefix + ".ins")

        # Run first iteration using xs
        self.driver.run_SHELXTL_command(cmd="xs")
        shutil.copy(os.path.join(self.path_to_ins, self.output_prefix + ".res"), os.path.join(self.path_to_ins, self.output_prefix + ".ins"))

        # Read in and run initial SHELXTL file
        ins_file = self.driver.get_ins_file()
        ins_file.remove_command('L.S.')
        ins_file.add_command('L.S.', [str(self.least_squares_iterations)])
        self.utils = OptimizerUtils(ins_file, self.bond_lengths, self.mixing_pairs, self.use_ml_model)
        self.history = OptimizerHistory(self.driver, self.utils, ins_file, self.score_weighting, self.max_n_leaves)

        # Optimization
        self.run_step(self.optimizer_steps.identify_sites)
        self.run_step(self.optimizer_steps.switch_elements)
        self.run_step(self.optimizer_steps.try_site_mixing)

        self.run_step(self.optimizer_steps.change_occupancy)
        self.run_step(self.optimizer_steps.try_exti)
        self.run_step(self.optimizer_steps.try_anisotropy)
        self.run_step(self.optimizer_steps.use_suggested_weights)
        self.run_step(self.optimizer_steps.use_suggested_weights)

        self.driver.run_SHELXTL(self.history.get_best_history()[-1].ins_file)
        print "Done with optimization"
        results_path = os.path.join(self.path_to_ins, "optimizer_results")
        if not os.path.exists(results_path):
            os.mkdir(results_path)
        with open(os.path.join(results_path, "report.txt"), 'w') as f:
            if self.history.head.r1 > 0.1:
                f.write("High initial R1 score, there may be something wrong with the site assigments (which are actually sites)\n")
            if self.history.get_best_history()[-1].r1 > 0.1:
                f.write("High final R1 score, the optimization may not have been successful\n")
            f.write(self.utils.get_report())
            print "Report on optimization process written to " + os.path.join(results_path, "report.txt")
            print "Output res files for top {} results saved in ".format(self.n_results) + results_path
        self.generate_graph(os.path.join(results_path, "optimization_graph"))
        print "Graph of optimization process saved to " + os.path.join(results_path, "optimization_graph.pdf")
        sorted_leaves = self.history.head.get_sorted_leaves()
        for i in range(0, min(self.n_results, len(self.history.leaves))):
            with open(os.path.join(results_path, "{}.res".format(min(self.n_results, len(self.history.leaves)))), 'w') as f:
                f.write(sorted_leaves[i-1].res_file.filetxt)
        # bonds = self.utils.get_bonds(self.driver, self.history.get_best_history()[-1].res_file)
        # bond_by_atom = defaultdict(lambda: [])
        # for bond in bonds:
        #     bond_by_atom[bond[0]].append((bond[1], bond[2]))
        #     bond_by_atom[bond[1]].append((bond[0], bond[2]))
        # from citrination_client import CitrinationClient
        # from pymatgen.core import Element
        # for atom, bonds in bond_by_atom.items():
        #     shortest = sorted(bonds, key=lambda tup: tup[1])[0]
        #     print atom, shortest
        #     candidate = {"Element 1": re.sub("\d", "", atom), "Element 2": re.sub("\d", "", shortest[0]), "formula": self.history.get_best_history()[-1].res_file.get_analytic_formula()}
        #     result = CitrinationClient(os.environ["CITRINATION_API_KEY"]).predict("680", candidate)["candidates"][0]["Bond length"]
        #     print "ml model:", shortest[1] - result[0], result
        #     print "sum of radii:", shortest[1] - (Element(re.sub("\d", "", atom)).atomic_radius + Element(re.sub("\d", "", shortest[0])).atomic_radius), Element(re.sub("\d", "", atom)).atomic_radius + Element(re.sub("\d", "", shortest[0])).atomic_radius
        #     print "#"*50

        # print map(lambda tup: (tup[0], sorted(tup[1])[0]), bond_by_atom.items())


    def run_step(self, step):
        """
        Run a single step in the optimization
        :param step: Which optimization step to try
        :return:
        """
        for leaf in self.history.leaves:
            step(leaf)

    def generate_graph(self, output_file):
        self.history.generate_graph(output_file)

    def check_inputs(self):
        if self.bond_lengths is not None:
            assert type(self.bond_lengths) == list, "bond_lengths argument should be a list of tuples of length 3"
            for bl in self.bond_lengths:
                assert type(bl) == tuple, "bond_lengths argument should be a list of tuples of length 3"
                assert len(bl) == 3, "bond_lengths argument should be a list of tuples of length 3"

        if self.mixing_pairs is not None:
            assert type(self.mixing_pairs) == list, "mixing_pairs argument should be a list of tuples of length 2"
            for mp in self.mixing_pairs:
                assert type(mp) == tuple, "mixing_pairs argument should be a list of tuples of length 2"
                assert len(mp) == 2, "mixing_pairs argument should be a list of tuples of length 2"

        try:
            float(self.r1_similarity_threshold)
        except ValueError:
            print("r1_similarity_threshold should be numerical")
        assert 0 <= self.r1_similarity_threshold <= 1.0, "r1_similarity_threshold should be a real number between 0.0 and 1.0"

        try:
            float(self.r1_threshold)
        except ValueError:
            print("r1_threshold should be numerical")
        assert 0 <= self.r1_threshold <= 1.0, "r1_threshold should be a real number between 0.0 and 1.0"

        try:
            float(self.occupancy_threshold)
        except ValueError:
            print("occupancy_threshold should be numerical")
        assert 0 <= self.occupancy_threshold <= 1.0, "occupancy_threshold should be a real number between 0.0 and 1.0"

        try:
            int(self.least_squares_iterations)
        except ValueError:
            print("least_squares_iterations should be numerical")
        assert 0 < self.least_squares_iterations, "least_squares_iterations should be a positive integer"
        assert int(self.least_squares_iterations) == self.least_squares_iterations, \
            "least_squares_iterations should be a positive integer"

        try:
            float(self.score_weighting)
        except ValueError:
            print("score_weighting should be numerical")
        assert 0 <= self.score_weighting <= 1.0, "score_weighting should be a real number between 0.0 and 1.0"

        try:
            int(self.max_n_leaves)
        except ValueError:
            print("max_n_leaves should be numerical")
        assert 0 < self.max_n_leaves, "max_n_leaves should be a positive integer"
        assert int(self.max_n_leaves) == self.max_n_leaves, "max_n_leaves should be a positive integer"

        try:
            int(self.n_results)
        except ValueError:
            print("n_results should be numerical")
        assert 0 < self.n_results, "n_results should be a positive integer"
        assert int(self.n_results) == self.n_results, "n_results should be a positive integer"



