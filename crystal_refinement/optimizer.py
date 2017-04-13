from SHELXDriver import SHELXDriver
import os, re, copy, math, itertools, random
import numpy as np
import shutil
from pymatgen.io.cif import CifParser
from pymatgen.core.composition import Element
import utils
from OptimizerHistory import OptimizerHistory


class Optimizer:
    """
    Class for performing single crystal refinement
    """
    def __init__(self, r1_similarity_threshold=0.0025, occupancy_threshold=0.02, r1_threshold=0.1,
                 least_squares_iterations=10, use_ml_model=False):
        """
        :param r1_similarity_threshold: If r1 scores for multiple options are within this similarity threshold, the
            optimizer will branch and explore all the options.
        :param occupancy_threshold: Minimum deviation in occupancy or site mixing to apply those steps
        :param r1_threshold: Threshold to denote a high r1 score, which triggers an alternate path in the optimizer
            where bond lengths drive the optimization instead of r1
        :param use_ml_model: Whether to use the bond length ml model to estimate bond lengths
        """
        self.r1_similarity_threshold = r1_similarity_threshold
        self.r1_threshold = r1_threshold
        self.occupancy_threshold = occupancy_threshold
        self.least_squares_iterations = least_squares_iterations
        self.use_ml_model = use_ml_model

    def run(self, path_to_xl, path_to_xs, ins_path, input_prefix, output_prefix, use_wine=False, annotate_graph=False):
        """
        Method to run the optimization
        :param path_to_xl: path to xl executable
        :param path_to_xs: path to xs executable
        :param ins_path: path to ins file output by xprep
        :param input_prefix: prefix of ins file (eg for file.ins, the prefix would be "file")
        :param output_prefix: prefix of the result file that the optimizer will output
        :return:
        """

        # Copy ins and hkl file to output prefix
        os.chdir(ins_path)
        shutil.copy(os.path.join(ins_path, input_prefix + ".hkl"), os.path.join(ins_path, output_prefix + ".hkl"))
        shutil.copy(os.path.join(ins_path, input_prefix + ".ins"), os.path.join(ins_path, output_prefix + ".ins"))
        self.annotate_graph = annotate_graph
        self.driver = SHELXDriver(ins_path=ins_path, prefix=output_prefix, path_to_xl=path_to_xl, path_to_xs=path_to_xs, use_wine=use_wine)

        # Check that the ins file is direct from xprep, without having been run before
        f = open(output_prefix + ".ins")
        assert len(f.readlines()) < 15, "Error: Must run optimizer directly on output from xprep, without other changes"

        # Run first iteration using xs
        self.driver.run_SHELXTL_command(cmd="xs")
        shutil.copy(os.path.join(ins_path, output_prefix + ".res"), os.path.join(ins_path, output_prefix + ".ins"))

        # Read in and run initial SHELXTL file
        ins_file = self.driver.get_ins_file()
        ins_file.remove_command('L.S.')
        ins_file.add_command('L.S.', [str(self.least_squares_iterations)])
        self.history = OptimizerHistory(self.driver, ins_file)

        # Optimization
        self.run_step(self.identify_sites)
        self.run_step(self.switch_elements)

        # There is currently an inherent assumption that switch elements will never be called after site mixing, since
        # after site mixing, the indices don't line up anymore
        self.run_step(self.try_site_mixing)
        self.run_step(self.change_occupancy)
        self.run_step(self.try_exti)
        self.run_step(self.try_anisotropy)
        self.run_step(self.use_suggested_weights)
        self.run_step(self.use_suggested_weights)

        self.driver.run_SHELXTL(self.history.get_best_history()[-1].ins_file)
        print "Done with optimization"

    def run_step(self, step):
        for leaf in self.history.leaves:
            step(leaf)

    def get_bonds(self, ins_file):
        """
        Get the bonds defined in the given ins_file

        :param ins_file: SHELXFile object
        :return res: List of bond tuples (element 1, element 2, bond length)
        """
        ins_file.add_command("ACTA")
        self.driver.run_SHELXTL(ins_file)
        ins_file.remove_command("ACTA")

        with open(self.driver.cif_file) as f:
            cif_file = CifParser.from_string(f.read())

        cif_dict = cif_file.as_dict().values()[0]
        return zip(cif_dict["_geom_bond_atom_site_label_1"], cif_dict["_geom_bond_atom_site_label_2"],
                   [float(x.replace("(", "").replace(")", "")) for x in cif_dict["_geom_bond_distance"]])

    def get_shortest_bond(self, ins_file):
        return sorted([utils.get_ideal_bond_length(el.capitalize(), el.capitalize()) for el in ins_file.elements])[0]

    def identify_sites(self, initial):
        if initial.r1 > self.r1_threshold:
            self.identify_sites_by_bond_length(initial)
        else:
            self.try_add_q(initial)
            for leaf in initial.get_leaves():
                self.try_remove_site(leaf)

    def identify_sites_by_bond_length(self, initial):
        """
        If the initial r1 is large, use bond lengths instead of the r1 score as the criteria for identifying which
        electron density peaks are atom sites
        """
        ins_file = initial.get_res()
        shortest_possible_bond = self.get_shortest_bond(ins_file)
        prev_iteration = initial
        for i in range(5):
            ins_file = prev_iteration.get_res()
            to_delete = set()
            for bond in sorted(self.get_bonds(ins_file), key=lambda tup: tup[2]):
                # Threshold on how short the bonds are
                if bond[2] / shortest_possible_bond < 0.5:
                    a1_num = int(re.search('\d+', bond[0]).group(0))
                    a2_num = int(re.search('\d+', bond[1]).group(0))
                    to_delete.add(max(a1_num, a2_num))
            if len(to_delete) == 0:
                break
            ins_file.remove_sites_by_number(to_delete)
            prev_iteration = self.history.run_and_save(ins_file, prev_iteration)

    def try_add_q(self, initial):
        """
        Try adding q peaks to main crystal sites if it decreases R value
        :return:
        """
        annotation = None
        if self.annotate_graph:
            annotation = "Added q peak"
        ins_file = initial.get_res()
        # This threshold could be scaled based on the potential atoms
        if ins_file.q_peaks[0].electron_density > 50:
            ins_file.move_q_to_crystal()
            # Find best element for new site
            num_elems = len(ins_file.elements)
            iterations = []
            for elem in range(1, num_elems + 1):
                ins_file.change_element(len(ins_file.crystal_sites) - 1, elem)
                iteration = self.history.run_iter(ins_file, initial, annotation)
                if iteration is not None:
                    iterations.append(iteration)
            iterations.sort(key=lambda i: i.r1)
            best_iter = iterations[0]
            displacements = [x.displacement for x in best_iter.crystal_sites]
            if best_iter.r1 < initial.r1 and displacements[-1] < (np.mean(displacements[:-1]) + 2.0 * np.std(displacements[:-1])):
                self.history.save([iterations[0]])
                for leaf in initial.get_leaves():
                    #   If adding one peak helped, recursively try adding another peak until it stops helping
                    self.try_add_q(leaf)

    def try_remove_site(self, initial):
        """
        Remove crystal sites if they result in bond distances that are too short
        :param ml_model: if True, a machine learning model is used to predict the correct bond length.
                         if False, the ideal bond length is approximated as the sum of the atomic radii
        :return:
        """

        ins_file = initial.get_res()
        r_penalty = 1.1
        bonds = self.get_bonds(ins_file)
        threshold = 0.1
        while True:
            ins_file = initial.get_res()
            to_delete = set()
            for a1, a2, distance in bonds:
                ideal_distance = utils.get_ideal_bond_length(a1, a2, self.use_ml_model)
                # if the distance is too small, remove the lower density atom
                if (ideal_distance - distance) / ideal_distance > threshold:
                    a1_num = int(re.search('\d+', a1).group(0))
                    a2_num = int(re.search('\d+', a2).group(0))
                    to_delete.add(max(a1_num, a2_num))
            # no sites removed
            if len(to_delete) == 0:
                break
            ins_file.remove_sites_by_number(to_delete)
            if self.annotate_graph:
                cur_iter = self.history.run_iter(ins_file, initial, "Removed {} site(s)".format(len(to_delete)))
            else:
                cur_iter = self.history.run_iter(ins_file, initial)
            if cur_iter is not None and cur_iter.r1 < initial.r1 * r_penalty:
                self.history.save(cur_iter)
                break
            threshold *= 1.1

    def switch_elements(self, initial):
        ins_file = initial.get_res()

        # Want to make changes from largest displacement to smallest
        sorted_sites = sorted(ins_file.crystal_sites, key=lambda s: -s.displacement)
        order = [s.site_number - 1 for s in sorted_sites]
        num_elems = len(ins_file.elements)
        # print "\n".join([" ".join(s.write_line()) for s in initial.ins_file.crystal_sites])
        for i in order:
            for prev_iter in initial.get_leaves():
                ins_file = prev_iter.get_res()
                iterations = []
                for elem in range(1, num_elems+1):
                    ins_file.change_element(i, elem)
                    if self.annotate_graph:
                        prev = prev_iter.res_file.crystal_sites[i].name
                        cur = ins_file.crystal_sites[i].name
                        iteration = self.history.run_iter(ins_file, prev_iter, "Changed {} to {}".format(prev, cur))
                    else:
                        iteration = self.history.run_iter(ins_file, prev_iter)
                    if iteration is not None:
                        iterations.append(iteration)
                iterations.sort(key=lambda i: i.r1)
                for iteration in iterations:
                    if iteration.r1 - iterations[0].r1 < self.r1_similarity_threshold:
                        self.history.save(iteration)

    def change_occupancy(self, initial):
        ins_file = initial.get_res()

        # Want to make changes from largest displacement to smallest
        displacements = map((lambda x: x.displacement), ins_file.crystal_sites)

        for i, displacement in sorted(enumerate(displacements), key=lambda tup: -tup[1]):
            # don't change occupancy of mixed sites
            if ins_file.crystal_sites[i].site_number in ins_file.mixed_site_numbers:
                continue
            for prev_iter in initial.get_leaves():
                ins_file = prev_iter.get_res()

                # only change occupancy if displacement is >= 2 std deviations away from mean
                mean = np.mean(displacements[:i] + displacements[i+1:])
                std = np.std(displacements[:i] + displacements[i+1:])
                if abs(displacement - mean) / std < 2.0:
                    break

                ins_file.add_variable_occupancy(i)
                if self.annotate_graph:
                    site = ins_file.crystal_sites[i].name
                    iteration = self.history.run_iter(ins_file, initial, "Added variable occupancy for {}".format(site))
                else:
                    iteration = self.history.run_iter(ins_file, initial)

                # If changing the occupancy decreased r1, decreased the displacement, and resulted in an occupancy
                # that satisfies the threshold, add it to the history
                if iteration is not None \
                        and iteration.r1 < prev_iter.r1 \
                        and float(iteration.res_file.fvar_vals[-1]) < (1 - self.occupancy_threshold) \
                        and iteration.res_file.crystal_sites[i].displacement < displacement:
                    self.history.save(iteration)

    def try_site_mixing(self, initial):
        """
        Try allowing site mixing, where a single crystal site might have mixed occupancy between two different elements
        Only handles pair-wise mixing.
        :param driver: SHELX driver
        """
        ins_file = initial.get_res()
        element_list = [Element(el.capitalize()) for el in ins_file.elements]
        pairs = []
        probability_threshold = 2E-4

        # For all elements in compound, calculate substitution probabilities
        # If substitution probability is > probability_threshold, then save it to pairs list.
        for i1, i2 in itertools.combinations(range(len(element_list)), 2):
            e1 = element_list[i1]
            e2 = element_list[i2]
            sp = utils.get_substitution_probability(e1, e2)
            if sp > probability_threshold:
                pairs.append(([i1, i2], sp))

        # Sort pairs by substitution probability (largest to smallest)
        pairs = [tup[0] for tup in sorted(pairs, key=lambda tup: -tup[1])]

        # Keep adding site mixing until we've already tried adding site mixing at the top priority sites
        tried = set()
        self.do_site_mixing(initial, tried, pairs)

    def do_site_mixing(self, initial, tried, pairs):
        ins_file = initial.get_res()
        bonds = self.get_bonds(ins_file)
        mixing_priority = utils.site_mixing_priority(bonds)
        # In case of ties, find all top tied priorities
        top_priority_score = mixing_priority[0][1]
        top_priority = [priority[0] for priority in mixing_priority if priority[1] - top_priority_score < 0.5]

        # For each of these top priorities, try site mixing
        if all([i in tried for i in top_priority]):
            return
        for i in top_priority:
            if i in tried:
                continue
            for prev_iter in initial.get_leaves():
                tried.add(i)
                iterations = []
                for pair in pairs:
                    ins_file = prev_iter.get_res()
                    ins_file.add_site_mixing(site_number=i, mixing_element_indices=pair)
                    if self.annotate_graph:
                        mix = "{} and {}".format(ins_file.elements[pair[0]], ins_file.elements[pair[1]])
                        iteration = self.history.run_iter(ins_file, prev_iter, "Mixing {} on site {}".format(mix, i))
                    else:
                        iteration = self.history.run_iter(ins_file, prev_iter)
                    if iteration is not None:
                        occupancy_var = float(iteration.res_file.fvar_vals[-1])
                        # Only include occupancies that are actually split
                        if occupancy_var > self.occupancy_threshold and occupancy_var < (1 - self.occupancy_threshold):
                            iterations.append(iteration)
                if len(iterations) == 0:
                    continue
                iterations.sort(key=lambda i: i.r1)
                best_r1 = min([prev_iter.r1, iterations[0].r1])
                if prev_iter.r1 - best_r1 < self.r1_similarity_threshold:
                    prev_iter.propagate()
                if iterations[0].r1 - best_r1 < self.r1_similarity_threshold:
                    self.history.save(iterations[0])
                for leaf in prev_iter.get_leaves():
                    self.do_site_mixing(leaf, tried.union(set(top_priority)), pairs)


    def try_anisotropy(self, initial):
        """
        Test if adding anisotropy reduces R1 value.  If it does, do so.
        :return:
        """

        ins_file = initial.get_res()

        #  Try with anisotropy
        ins_file.add_anisotropy()
        if self.annotate_graph:
            iteration = self.history.run_iter(ins_file, initial, "Added anisotropy")
        else:
            iteration = self.history.run_iter(ins_file, initial)
        #  If anisotropy helped, add it to the history
        if iteration is not None and iteration.r1 < initial.r1:
            self.history.save(iteration)

    def try_exti(self, initial):
        """
        Test if adding extinguishing reduces R1 value.  If it does, do so.
        :return:
        """
        ins_file = initial.get_res()

        #  Try with extinguishing
        ins_file.add_exti()
        if self.annotate_graph:
            iteration = self.history.run_iter(ins_file, initial, "Added extinction")
        else:
            iteration = self.history.run_iter(ins_file, initial)

        #  If exti helped, add it to the history
        if iteration is not None and iteration.r1 < initial.r1:
            self.history.save(iteration)

    def use_suggested_weights(self, initial):
        """
        Stop re-initializing weights each time--use previously suggested weights
        :return:
        """
        ins_file = initial.get_res()

        #  Try with extinguishing
        ins_file.remove_command("WGHT")
        ins_file.commands.append(("WGHT", ins_file.suggested_weight_vals))
        if self.annotate_graph:
            iteration = self.history.run_iter(ins_file, initial, "Used suggested weights")
        else:
            iteration = self.history.run_iter(ins_file, initial)

        #  If exti helped, add it to the history
        if iteration is not None and iteration.r1 < initial.r1:
            self.history.save(iteration)


########################################################################################################################


def test_all(path_to_SXTL_dir, ins_folder, input_prefix="absfac1", output_prefix="temp", use_wine=False,
             generate_graph=False, annotate_graph=False, graph_path=""):
    subdirs = os.listdir(ins_folder)
    for dirname in subdirs:
        if dirname[0] != ".":
            print dirname
            test_single(path_to_SXTL_dir, os.path.join(ins_folder, dirname), input_prefix, output_prefix, use_wine,
                        generate_graph, annotate_graph, graph_path)


def test_single(path_to_SXTL_dir, dirname, input_prefix="absfac1", output_prefix="temp", print_files=False, use_wine=False,
                generate_graph=False, annotate_graph=False, graph_path=""):
    try:
        ins_path = os.path.join(dirname, "work") + "/"
        for filename in os.listdir(os.path.join(dirname, "Anton")):
            if ".hkl" in filename:
                shutil.copy(os.path.join(dirname, "Anton", filename),
                            os.path.join(ins_path, input_prefix + ".hkl"))
            if ".res" in filename:
                final_res = os.path.join(dirname, "Anton", filename)
    except Exception:
        try:
            ins_path = dirname + "/"
            open(os.path.join(dirname, "1.hkl"))
            open(os.path.join(dirname, "1.ins"))
            final_res = os.path.join(dirname, "result.res")
            open(final_res)
        except Exception, e:
            print "File structure failure", e
            print "~" * 50
            return


    opt = run_single(path_to_SXTL_dir, ins_path, input_prefix, output_prefix, use_wine,
                     generate_graph, annotate_graph, graph_path)

    best_history = opt.history.get_best_history()

    r1_tol = 2e-4
    anton_r1 = float(re.search("REM R1 =  (\d\.\d+)", open(final_res).read()).group(1))
    print len(opt.history.leaves), "path(s) tried"
    print "Initial r1 = {}".format(best_history[0].r1)
    print "Optimizer r1 = {}".format(best_history[-1].r1)
    print "Reference r1 = {}".format(anton_r1)
    if best_history[-1].r1 - anton_r1 > r1_tol:
        print "Not success!"
    if print_files:
        print "Optimizer final result:"
        print open(os.path.join(ins_path, output_prefix + ".res")).read()
        print "Reference final result:"
        print open(final_res).read()
    print "~" * 50


def run_single(path_to_SXTL_dir, ins_path, input_prefix="absfac1", output_prefix="temp", use_wine=False,
               generate_graph=False, annotate_graph=False, graph_path=""):
    opt = Optimizer()
    opt.run(os.path.join(path_to_SXTL_dir, "xl.exe"), os.path.join(path_to_SXTL_dir, "xs.exe"), ins_path, input_prefix, output_prefix, use_wine=use_wine,
            annotate_graph=annotate_graph)
    if generate_graph:
        print os.path.join(graph_path, os.path.basename(ins_path))
        opt.history.head.generate_graph(os.path.join(graph_path, os.path.basename(ins_path)))
    return opt


def run_all(path_to_SXTL_dir, ins_folder, input_prefix="absfac1", output_prefix="temp", use_wine=False,
            generate_graph=False, annotate_graph=False, graph_path=""):
    subdirs = os.listdir(ins_folder)
    print subdirs
    for dirname in subdirs:
        if dirname[0] != ".":
            print dirname
            opt = run_single(path_to_SXTL_dir, os.path.join(ins_folder, dirname), input_prefix, output_prefix, use_wine,
                             generate_graph, annotate_graph, graph_path)
            best_history = opt.history.get_best_history()
            print len(opt.history.leaves), "path(s) tried"
            print "Initial r1: {}".format(best_history[0].r1)
            print "Final r1: {}".format(best_history[-1].r1)

def main():
    path_to_SXTL_dir = "/Users/eantono/Documents/program_files/xtal_refinement/SXTL/"
    # ins_folder = "/Users/eantono/Documents/project_files/xtal_refinement/4-2-1-4 INS and HKL files"
    ins_folder = "/Users/eantono/Documents/project_files/xtal_refinement/!UNSEEN 4-2-1-4/"
    subdir = "!solid_solution-Nd4Mn2CdSi2.5Ge1.5"
    graph_output_path = "/Users/eantono/Documents/src/xtal_refinement/output"
    # path_to_SXTL_dir = "/Users/julialing/Documents/GitHub/crystal_refinement/shelxtl/SXTL/"
    # ins_folder = "/Users/julialing/Documents/DataScience/crystal_refinement/single_crystal_data/"

    # test_all(path_to_SXTL_dir, ins_folder, input_prefix="1", use_wine=True,
    #   generate_graph=True, annotate_graph=True, graph_path=graph_output_path)
    # test_single(path_to_SXTL_dir, os.path.join(ins_folder, subdir), "1", print_files=True, use_wine=True,
    #   generate_graph=True, annotate_graph=True, graph_path=graph_output_path)
    run_all(path_to_SXTL_dir, ins_folder, input_prefix="1", use_wine=True,
      generate_graph=True, annotate_graph=True, graph_path=graph_output_path)
    # run_single(path_to_SXTL_dir, os.path.join(ins_path, "work/"), "absfac1",
    #   generate_graph=True, annotate_graph=True, graph_path=graph_output_path)

if __name__ == "__main__":
    main()


