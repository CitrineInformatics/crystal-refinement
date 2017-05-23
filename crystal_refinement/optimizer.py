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
                 write_graph=True, bond_lengths=None, mixing_pairs=None, use_ml_model=False,
                 r1_similarity_threshold=0.0075, occupancy_threshold=0.02, r1_threshold=0.1, score_weighting=0.8,
                 max_n_leaves=50, least_squares_iterations=4, n_results=10):
        """
        :param path_to_xl: path to xl executable (including executable file name)
        :param path_to_xs: path to xs executable (including executable file name)
        :param path_to_ins: path to directory containing *.ins file output by xprep (not including *.ins file name)
        :param input_prefix: prefix of ins file (e.g. for file.ins, the prefix would be "file")
        :param output_prefix: prefix of the result file that the optimizer will output
        :param use_wine: In order to run Windows executables on a mac, the wine app can be used.
        :param write_graph: If True, will write tree graph of optimizer paths to an image file
        :param bond_lengths: A list of prescribed ideal bond lengths in Angstroms.  Should be in the form of triples.
            (e.g. [(Ag, Ag, 2.85)] ).  If a bond length is not specified, then it is calculated via either a machine
            learning model or based on the covalent radii.
        :param mixing_pairs: A list of acceptable elements to allow to co-occupy a site.  If no such list is provided,
            then acceptable mixing pairs are determined based on pymatgen information.
        :param use_ml_model: Whether to use the bond length machine learning model to estimate bond lengths.
            If True, a machine learning model hosted at citrination.com is queried to determine the bond lengths.
                Using this option requires obtaining a free Citrination user account and setting the CITRINATION_API_KEY
                environment variable with your account API key.
                If this machine learning model has high uncertainty for the bond lengths queried, the optimizer defaults
                back to the covalent radii model.
            If False, a simpler bond length model based on covalent radii is used.
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
        self.write_graph = write_graph
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
        self.history = OptimizerHistory(self.driver, ins_file, self.score_weighting, self.max_n_leaves)
        self.utils = OptimizerUtils(ins_file, self.bond_lengths, self.mixing_pairs, self.use_ml_model)

        # Optimization
        self.run_step(self.optimizer_steps.identify_sites)
        self.run_step(self.optimizer_steps.switch_elements)
        # self.run_step(self.switch_elements)
        # There is currently an inherent assumption that switch elements will never be called after site mixing, since
        # after site mixing, the indices don't line up anymore
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
        for i in range(1, min(self.n_results, len(self.history.leaves)) + 1):
            with open(os.path.join(results_path, "{}.res".format(i)), 'w') as f:
                f.write(self.history.get_best_history()[-1*i].res_file.filetxt)
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

########################################################################################################################


def test_all(path_to_SXTL_dir, ins_folder, input_prefix="absfac1", output_prefix="temp", use_wine=False, print_files=False,
             generate_graph=False, truncated_graph=False, graph_path=""):
    subdirs = os.listdir(ins_folder)
    for dirname in subdirs:
        if dirname[0] != ".":
            print dirname
            test_single(path_to_SXTL_dir, os.path.join(ins_folder, dirname), input_prefix, output_prefix, use_wine,
                        print_files, generate_graph, truncated_graph, graph_path)


def test_single(path_to_SXTL_dir, dirname, input_prefix="absfac1", output_prefix="temp", use_wine=False, print_files=False,
                generate_graph=False, truncated_graph=False, graph_path="", result_filename=None):
    if "INS-HKL" in os.listdir(dirname):
        final_res = ""
        ins_path = os.path.join(dirname, "INS-HKL")
        if result_filename is not None:
            final_res = os.path.join(ins_path, result_filename)
            input_prefix = input_prefix
        else:
            for filename in sorted(os.listdir(ins_path), key=lambda name: os.path.getctime(os.path.join(ins_path, name))):
                if ".res" in filename:
                    final_res = os.path.join(ins_path, filename)
                    break
            input_prefix = os.path.basename(final_res).split(".")[0]
        graph_name = os.path.basename(dirname)
        ins_from_result(ins_path, result_file=os.path.basename(final_res), input_prefix=input_prefix)
    else:
        try:
            ins_path = os.path.join(dirname, "work")
            graph_name = os.path.basename(dirname)
            for filename in os.listdir(os.path.join(dirname, "Anton")):
                if ".hkl" in filename:
                    shutil.copy(os.path.join(dirname, "Anton", filename),
                                os.path.join(ins_path, input_prefix + ".hkl"))
                if ".res" in filename:
                    final_res = os.path.join(dirname, "Anton", filename)
        except Exception:
            try:
                ins_path = dirname
                graph_name = os.path.basename(dirname)
                open(os.path.join(dirname, "1.hkl"))
                open(os.path.join(dirname, "1.ins"))
                final_res = os.path.join(dirname, "result.res")
                open(final_res)
            except Exception, e:
                print "File structure failure", e
                print "~" * 50
                return

    # try:

    start = time.time()
    opt = run_single(path_to_SXTL_dir, ins_path, input_prefix, output_prefix, use_wine,
                     generate_graph, truncated_graph, graph_path, graph_name)
    runtime = time.time() - start
    # except Exception, e:
    #     print "Optimizer failure:", e
    #     print "~" * 50
    #     return

    best_history = opt.history.get_best_history()

    r1_tol = 2e-4
    anton_r1 = float(re.search("REM R1 =  (\d\.\d+)", open(final_res).read()).group(1))
    print len(opt.history.leaves), "path(s) tried"
    print "Initial r1 = {}".format(best_history[0].r1)
    print "Optimizer r1 = {}".format(best_history[-1].r1)
    print "Reference r1 = {}".format(anton_r1)
    print "Runtime = {}s".format(runtime)
    if best_history[-1].r1 - anton_r1 > r1_tol:
        print "Not success!"
    if print_files:
        print "Optimizer final result:"
        print open(os.path.join(ins_path, output_prefix + ".res")).read()
        print "Reference final result:"
        print open(final_res).read()
    print "~" * 50


def run_single(path_to_SXTL_dir, ins_path, input_prefix="absfac1", output_prefix="temp", use_wine=False,
               generate_graph=False, truncated_graph=False, graph_path="", graph_name="out"):
    opt = Optimizer(os.path.join(path_to_SXTL_dir, "xl.exe"), os.path.join(path_to_SXTL_dir, "xs.exe"), ins_path, input_prefix,
            output_prefix, use_wine=use_wine, n_results=0, use_ml_model=True)
    opt.run()
    if generate_graph:
        if truncated_graph:
            opt.history.head.update_dead_branches()
            opt.history.head.generate_truncated_graph(os.path.join(graph_path, graph_name + "_trunc"))
        else:
            opt.history.head.generate_graph(os.path.join(graph_path, graph_name))
    return opt


def run_all(path_to_SXTL_dir, ins_folder, input_prefix="absfac1", output_prefix="temp", use_wine=False,
            generate_graph=False, truncated_graph=False, graph_path=""):
    subdirs = os.listdir(ins_folder)
    print subdirs
    for dirname in subdirs:
        if dirname[0] != ".":
            opt = run_single(path_to_SXTL_dir, os.path.join(ins_folder, dirname), input_prefix, output_prefix, use_wine,
                             generate_graph, truncated_graph, graph_path)
            best_history = opt.history.get_best_history()
            print len(opt.history.leaves), "path(s) tried"
            print "Initial r1: {}".format(best_history[0].r1)
            print "Final r1: {}".format(best_history[-1].r1)


def ins_from_result(folder_path, result_file="result.res", input_prefix="1"):
    from collections import defaultdict
    lines = defaultdict(lambda: [])
    lines["TREF"].append("\n")


    with open(os.path.join(folder_path, result_file)) as f:
        for line in f:
            split = line.split()
            if len(split) == 1:
                lines[split[0]].append("\n")
            elif len(split) > 1:
                lines[split[0]].append(line[len(split[0]):])


    keys = ["TITL", "CELL", "ZERR", "LATT", "SYMM", "SFAC", "UNIT", "TEMP", "SIZE", "TREF", "HKLF", "END"]
    with open(os.path.join(folder_path, input_prefix + ".ins"), 'w') as f:
        for key in keys:
            values = lines[key]
            for val in values:
                f.write("{} {}".format(key, val))



def main():
    path_to_SXTL_dir = "/Users/eantono/Documents/program_files/xtal_refinement/SXTL/"
    # ins_folder = "/Users/eantono/Documents/project_files/xtal_refinement/4-2-1-4 INS and HKL files"
    # ins_folder = "/Users/eantono/Documents/project_files/xtal_refinement/!UNSEEN 4-2-1-4/"
    # ins_folder = "/Users/eantono/Documents/project_files/xtal_refinement/Organized_data1/EASY"
    # ins_folder = "/Users/eantono/Documents/project_files/xtal_refinement/Organized_data1/MIXING"
    ins_folder = "/Users/eantono/Documents/project_files/xtal_refinement/Organized_data2/EASY"
    # ins_folder = "/Users/eantono/Documents/project_files/xtal_refinement/Organized_data2/Partial occupancy"
    # ins_folder = "/Users/eantono/Documents/project_files/xtal_refinement/Organized_data2/mixing"
    subdir = "Ce4Mn2InGe4_mar1247"
    graph_output_path = "/Users/eantono/Documents/src/xtal_refinement/output"
    # path_to_SXTL_dir = "/Users/julialing/Documents/GitHub/crystal_refinement/shelxtl/SXTL/"
    # ins_folder = "/Users/julialing/Documents/DataScience/crystal_refinement/single_crystal_data/"
    # from shutil import copyfile
    # for subdir in os.listdir(ins_folder):
    #     if subdir[0] != ".":
    #         folder = os.path.join(ins_folder, subdir, "INS-HKL")
    #         copyfile(os.path.join(folder, "1.hkl"), os.path.join(folder, "result.hkl"))
    # quit()
    test_all(path_to_SXTL_dir, ins_folder, input_prefix="1", use_wine=True, print_files=False,
      generate_graph=True, truncated_graph=True, graph_path=graph_output_path)
    # test_single(path_to_SXTL_dir, os.path.join(ins_folder, subdir), "1", use_wine=True, print_files=True,
    #   generate_graph=True, truncated_graph = True, graph_path=graph_output_path, result_filename=None)

    # run_all(path_to_SXTL_dir, ins_folder, input_prefix="1", use_wine=True,
    #   generate_graph=True, graph_path=graph_output_path)
    # run_single(path_to_SXTL_dir, os.path.join(ins_folder, subdir), "1", use_wine=True,
    #   generate_graph=True, graph_path=graph_output_path)

if __name__ == "__main__":
    main()


