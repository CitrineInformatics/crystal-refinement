from __future__ import absolute_import
import os, shutil, re, ast

import crystal_refinement.OptimizerSteps as OptimizerSteps
from crystal_refinement.SHELX.SHELXDriver import SHELXDriver
from crystal_refinement.history.OptimizerHistory import OptimizerHistory
from crystal_refinement.utils.OptimizerCache import OptimizerCache
from argparse import ArgumentParser

try:
    import configparser
except ImportError:
    import ConfigParser as configparser


class Optimizer:
    """
    Class for performing single crystal refinement
    """

    def __init__(self,
                 input_prefix,
                 output_prefix,
                 input_directory,
                 path_to_xl,
                 path_to_xs,
                 use_wine=False,
                 bond_lengths=None,
                 mixing_pairs=None,
                 use_ml_model=False,
                 citrination_api_key=None,
                 ensure_identified_elements=True,
                 r1_similarity_threshold=0.0075,
                 occupancy_threshold=0.02,
                 r1_threshold=0.1,
                 overall_score_ratio_threshold=1.3,
                 score_weighting=0.5,
                 max_n_leaves=50,
                 n_results=10,
                 suppress_output=True,
                 log_output=False):
        """
        :param path_to_xl: path to xl executable (including executable file name)
        :param path_to_xs: path to xs executable (including executable file name)
        :param input_directory: path to directory containing *.ins file output by xprep (not including *.ins file name)
        :param input_prefix: prefix of ins file (e.g. for file.ins, the prefix would be "file")
        :param output_prefix: prefix of the result file that the optimizer will output
        :param use_wine: In order to run Windows executables on a mac, the wine app can be used.
        :param bond_lengths: A list of prescribed ideal bond lengths in Angstroms.  Should be in the form of triples.
            (e.g. [('Ag', 'Ag', 2.85)] ).  If a bond length is not specified, then it is calculated via either a machine
            learning model or based on the atomic radii.
        :param mixing_pairs: A list of acceptable elements to allow to co-occupy a site. (e.g. [('Ge', 'Ru')])
            If no such list is provided, then acceptable mixing pairs are determined based on pymatgen information.
        :param use_ml_model: Whether to use the bond length machine learning model to estimate bond lengths.
            If True, a machine learning model hosted at https://citrination.com is queried to determine the bond lengths.
                If this machine learning model has high uncertainty for the bond lengths queried, the optimizer defaults
                back to the atomic radii model.
            If False, a simpler bond length model based on atomic radii is used.
        :param citrination_api_key: API key in order to use the Citrination bond length model. This can be found at
            https://citrination.com/users/edit
        :param ensure_identified_elements: Whether to penalize the score a result if any elements identified in the
            initial .ins file are not present.
        :param r1_similarity_threshold: Difference threshold for r1 scores to allow for branching.  (Should be double
            in range (0.0, 1.0)).
        :param occupancy_threshold: Minimum deviation in occupancy or site mixing to trigger partial occupancy
            (Should be double in range (0.0, 1.0)).
        :param r1_threshold: Threshold to denote a high r1 score, which triggers an alternate path in the optimizer
            where bond lengths drive the assignment of sites instead of the r1 score (Should be a double in range (0.0, 1.0)).
        :param overall_score_ratio_threshold: Ratio threshold for overall score similarity to allow for branching
        :param score_weighting: Weighting of R1 score versus bond length score when calculating the overall score.
            A value of 1.0 corresponds to R1 only.  (Should be a double in range (0.0, 1.0)).
        :param max_n_leaves: Maximum number of optimization paths to follow at any given optimization step.
            A higher value may give a more optimal path, but will be more computationally expensive.
            (Should be positive integer).
        :param n_results: Number of final .res files to save.  For n_results=1, only the results file with the
            best final score (based on R1 and bond lengths) will be saved.  (Should be a non-negative integer).
        :param suppress_output: Whether to suppress the output of the SHELXTL commands
        :param log_output: Whether to print the logs during the optimization process
        """

        # Initialize parameters based on arguments
        self.path_to_xl = path_to_xl
        self.path_to_xs = path_to_xs
        self.input_directory = os.path.abspath(input_directory)
        self.input_prefix = input_prefix
        self.output_prefix = output_prefix
        self.use_wine = use_wine
        self.bond_lengths = bond_lengths
        self.mixing_pairs = mixing_pairs
        self.use_ml_model = use_ml_model
        self.citrination_api_key = citrination_api_key
        self.ensure_identified_elements = ensure_identified_elements
        self.overall_score_ratio_threshold = overall_score_ratio_threshold
        self.r1_similarity_threshold = r1_similarity_threshold
        self.r1_threshold = r1_threshold
        self.occupancy_threshold = occupancy_threshold
        self.score_weighting = score_weighting
        self.max_n_leaves = max_n_leaves
        self.n_results = n_results
        self.log_output = log_output

        self.validate_inputs()

        # Initialize objects to be used during run()
        self.driver = SHELXDriver(directory=self.input_directory, prefix=self.output_prefix, path_to_xl=self.path_to_xl,
                                  path_to_xs=self.path_to_xs, use_wine=self.use_wine, suppress_ouput=suppress_output)
        self.history = None
        self.cache = None

    def run(self):
        """
        Method to run the optimization

        :return:
        """

        # Copy ins and hkl file to output prefix
        # os.chdir(self.input_directory)

        shutil.copy(os.path.join(self.input_directory, self.input_prefix + ".hkl"),
                    os.path.join(self.input_directory, self.output_prefix + ".hkl"))
        shutil.copy(os.path.join(self.input_directory, self.input_prefix + ".ins"),
                    os.path.join(self.input_directory, self.output_prefix + ".ins"))

        # Check that the ins file is direct from xprep, without having been run before
        # f = open(self.output_prefix + ".ins")

        # Run first iteration using xs
        self.driver.run_SHELXTL_command(cmd="xs")
        shutil.copy(os.path.join(self.input_directory, self.output_prefix + ".res"),
                    os.path.join(self.input_directory, self.output_prefix + ".ins"))

        # Read in and run initial SHELXTL file
        ins_file = self.driver.get_ins_file()
        if self.log_output:
            print("Initial INS file:")
            print(ins_file.to_string())
        ins_file.remove_command('L.S.')
        ins_file.add_command('L.S.', ["4"])
        self.cache = OptimizerCache(ins_file,
                                    self.bond_lengths,
                                    self.mixing_pairs,
                                    self.use_ml_model,
                                    self.citrination_api_key)

        self.history = OptimizerHistory(self.driver, self.cache, ins_file, self.score_weighting, self.max_n_leaves)

        # Optimization
        self.run_step(OptimizerSteps.identify_sites)
        self.run_step(OptimizerSteps.switch_elements)
        self.history.clean_history()
        self.run_step(OptimizerSteps.change_occupancy)
        self.run_step(OptimizerSteps.try_exti)
        self.run_step(OptimizerSteps.try_anisotropy)
        pre_weight_leaves = self.history.get_leaves()
        self.run_step(OptimizerSteps.use_suggested_weights)
        self.run_step(OptimizerSteps.use_suggested_weights)

        for pre_weight_leaf in pre_weight_leaves:
            self.history.clean_history(1, pre_weight_leaf)

        self.run_step(OptimizerSteps.try_site_mixing)
        self.history.clean_history(criteria=["overall_score", "r1_only"])

        self.driver.run_SHELXTL(self.history.get_best_history()[-1].ins_file)
        print("Done with optimization")
        results_path = os.path.join(self.input_directory, "optimizer_results")
        if not os.path.exists(results_path):
            os.mkdir(results_path)
        with open(os.path.join(results_path, "report.txt"), 'w') as f:
            if self.history.head.r1 > 0.1:
                f.write(
                    "WARN: High initial R1 score, there may be something wrong with the assignment of sites to electron density peaks\n")
            if self.history.get_best_history()[-1].r1 > 0.1:
                f.write("WARN: High final R1 score, the optimization may not have been successful\n")
            f.write(self.cache.get_report())
            print("Report on optimization process written to " + os.path.join(results_path, "report.txt"))
            print("Output res files for top {} results saved in ".format(self.n_results) + results_path)
        self.generate_graph(os.path.join(results_path, "optimization_graph"))
        print("Graph of optimization process saved to " + os.path.join(results_path, "optimization_graph.pdf"))
        sorted_leaves = self.history.head.get_sorted_leaves()
        for i in range(0, min(self.n_results, len(self.history.get_leaves()))):
            with open(os.path.join(results_path, "{}.res".format(i)), 'w') as f:
                f.write(sorted_leaves[i].res_file.filetxt)

    def run_step(self, step):
        """
        Run a single step in the optimization
        :param step: Which optimization step to try
        :return:
        """
        for leaf in self.history.get_leaves():
            step(leaf, self)

    def generate_graph(self, output_file):
        """
        Generate a graphviz graph of the optimizer history.
        :param output_file: File location to output the graph
        :return:
        """
        self.history.generate_graph(output_file)

    def validate_inputs(self):
        """
        Make sure that the optimizer arguments are valid. This will throw an assertion error for any invalid argument.
        :return:
        """

        if self.use_ml_model:
            assert self.citrination_api_key is not None, "To use the machine learning bond length model, you must " \
                                                         "specify your Citrination API key, which is available at " \
                                                         "https://citrination.com/users/edit. If you don't have a " \
                                                         "Citrination account, you can create a free account at " \
                                                         "https://citrination.com"

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


if __name__ == "__main__":
    """
    Command line interface for the optimizer
    """
    arg_parser = ArgumentParser()

    # arg_parser.add_argument("-d", "--directory",
    #                         help="Path to the directory containing the initial ins file, without filename, "
    #                              "e.g. /path/to/file/", dest="directory")
    # arg_parser.add_argument("-i", "--input", dest="input_prefix",
    #                         help="Prefix for input file. e.g. if input file is input.ins, then the prefix is 'input'")
    # arg_parser.add_argument("-o", "--output", dest="output_prefix",
    #                         help="Prefix for output file. This is where result of optimization will be written")
    # arg_parser.add_argument("-c", "--config", dest="config_file",
    #                         help="Configuration file specifying additional parameters for the optimizer")

    arg_parser.add_argument("directory", help="Path to the directory containing the initial ins file, without filename,"
                                              " e.g. /path/to/file/")
    arg_parser.add_argument("input_prefix", help="Prefix for input file. e.g. if input file is input.ins, then the "
                                                 "prefix is 'input'")
    arg_parser.add_argument("output_prefix", help="Prefix for output file. This is where result of optimization will "
                                                  "be written")
    arg_parser.add_argument("config_file", help="Configuration file specifying additional parameters for the optimizer")

    args = arg_parser.parse_args()
    input_directory = args.directory
    input_prefix = args.input_prefix
    output_prefix = args.output_prefix

    config_parser = configparser.ConfigParser()
    config_parser.read(args.config_file)

    # SHELX Config section
    xl_path = config_parser.get("SHELX Config", "path_to_xl")
    xs_path = config_parser.get("SHELX Config", "path_to_xs")
    use_wine = config_parser.getboolean("SHELX Config", "use_wine", fallback=False)

    # Optimizer Config section
    r1_threshold = config_parser.getfloat("Optimizer config", "r1_threshold", fallback=0.1)
    r1_similarity_threshold = config_parser.getfloat("Optimizer config", "r1_similarity_threshold", fallback=0.0075)
    overall_score_ratio_threshold = config_parser.getfloat("Optimizer config", "overall_score_ratio_threshold",
                                                           fallback=1.3)
    occupancy_threshold = config_parser.getfloat("Optimizer config", "occupancy_threshold", fallback=0.02)
    max_n_leaves = config_parser.getfloat("Optimizer config", "max_n_leaves", fallback=50)

    # Bond Length Config section
    bond_lengths = config_parser.get("Bond Length Config", "bond_lengths", fallback=None)
    if bond_lengths is not None:
        bond_length_assert_message = "{} is an invalid configuration for bond_lengths. Please enter a list of tuples " \
                                     "in the form [('Ge', 'Mn', 2.5)]".format(bond_lengths)
        try:
            bond_lengths = ast.literal_eval(bond_lengths)
            assert type(bond_lengths) == list, bond_length_assert_message
            for bond in bond_lengths:
                assert type(bond) == tuple \
                       and len(bond) == 3 \
                       and type(bond[0]) == str \
                       and type(bond[1]) == str \
                       and type(bond[2]) == float, bond_length_assert_message

        except ValueError:
            assert False, bond_length_assert_message
    use_ml_model = config_parser.getboolean("Bond Length Config", "use_ml_model", fallback=False)
    citrination_api_key = config_parser.get("Bond Length Config", "citrination_api_key", fallback=None)

    # Mixing Pairs Config section
    mixing_pairs = config_parser.get("Mixing Pairs Config", "mixing_pairs", fallback=None)
    if mixing_pairs is not None:
        mixing_pairs_assert_message = "{} is an invalid configuration for bond_lengths. Please enter a list of tuples " \
                                      "in the form [('Au', 'Os')]".format(bond_lengths)
        try:
            mixing_pairs = ast.literal_eval(mixing_pairs)
            assert type(mixing_pairs) == list, mixing_pairs_assert_message
            for mixing_pair in mixing_pairs:
                assert type(mixing_pair) == tuple \
                       and len(mixing_pair) == 2 \
                       and type(mixing_pair[0]) == str \
                       and type(mixing_pair[1]) == str, mixing_pairs_assert_message

        except ValueError:
            assert False, mixing_pairs_assert_message

    # Score Config section
    score_weighting = config_parser.getfloat("Score Config", "score_weighting", fallback=1.0)
    ensure_identified_elements = config_parser.getboolean("Score Config", "ensure_identified_elements", fallback=True)

    # Logging Config section
    n_results = config_parser.getint("Logging Config", "n_results", fallback=10)
    suppress_output = config_parser.getboolean("Logging Config", "suppress_output", fallback=True)
    log_output = config_parser.getboolean("Logging Config", "log_output", fallback=False)

    # Create optimizer object and run on example
    opt = Optimizer(input_prefix=input_prefix,
                    output_prefix=output_prefix,
                    input_directory=input_directory,
                    path_to_xl=xl_path,
                    path_to_xs=xs_path,
                    use_wine=use_wine,
                    bond_lengths=bond_lengths,
                    mixing_pairs=mixing_pairs,
                    use_ml_model=use_ml_model,
                    citrination_api_key=citrination_api_key,
                    ensure_identified_elements=ensure_identified_elements,
                    r1_similarity_threshold=r1_similarity_threshold,
                    occupancy_threshold=occupancy_threshold,
                    r1_threshold=r1_threshold,
                    overall_score_ratio_threshold=overall_score_ratio_threshold,
                    score_weighting=score_weighting,
                    max_n_leaves=max_n_leaves,
                    n_results=n_results,
                    suppress_output=suppress_output,
                    log_output=log_output)
    opt.run()
