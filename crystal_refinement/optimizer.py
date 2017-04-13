from SHELXDriver import SHELXDriver
import os, re, copy, math, itertools, random
import numpy as np
import shutil
from pymatgen.io.cif import CifParser
from pymatgen.core.composition import Element
import utils


class Optimizer:
    """
    Class for performing single crystal refinement
    """
    def __init__(self):
        self.ins_history = []  # History of all previous ins files
        self.r1_history = []  # History of all previous R1 values

    def run(self, path_to_xl, path_to_xs, ins_path, input_prefix, output_prefix, use_wine=False):
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

        self.driver = SHELXDriver(ins_path=ins_path, prefix=output_prefix, path_to_xl=path_to_xl, path_to_xs=path_to_xs, use_wine=use_wine)

        # Check that the ins file is direct from xprep, without having been run before
        f = open(output_prefix + ".ins")
        assert len(f.readlines()) < 15, "Error: Must run optimizer directly on output from xprep, without other changes"

        # Run first iteration using xs
        self.driver.run_SHELXTL_command(cmd="xs")
        shutil.copy(os.path.join(ins_path, output_prefix + ".res"), os.path.join(ins_path, output_prefix + ".ins"))

        # Read in and run initial SHELXTL file
        ins_file = self.driver.get_ins_file()
        for cmd in ins_file.commands:
            if cmd[0] == "L.S.":
                ins_file.commands.remove(cmd)
                ins_file.commands.insert(0, ('L.S.', [str(10)]))
        self.run_iter(ins_file)

        # Optimization
        if self.r1_history[-1] > 0.1:
            self.identify_sites()
        else:
            self.try_add_q()
            self.try_remove_site()

        ins_file = self.driver.get_res_file()
        self.run_iter(ins_file)

        self.switch_elements()
        # There is currently an inherent assumption that switch elements will never be called after site mixing, since
        # after site mixing, the indices don't line up anymore

        self.try_site_mixing()


        self.change_occupancy()
        ins_file = self.driver.get_res_file()
        self.run_iter(ins_file)
        self.try_exti()
        self.try_anisotropy()
        self.use_suggested_weights()
        self.use_suggested_weights()

        print "Done with optimization"



    def run_iter(self, ins_file, ins_history=None, r1_history=None):
        """
        Run the given ins file through SHELXTL and record the file and resulting r1

        :param ins_file: SHELXFile object
        :return res: SHELXFile object
        """
        if ins_history is None:
            ins_history = self.ins_history

        if r1_history is None:
            r1_history = self.r1_history
        ins_history.append(copy.deepcopy(ins_file))
        res = self.driver.run_SHELXTL(ins_file)
        if res is not None:
            r1_history.append(res.r1)
        else:
            r1_history.append(np.infty)
        return res

    def get_bonds(self, ins_file):
        """
        Get the bonds defined in the given ins_file

        :param ins_file: SHELXFile object
        :return res: List of bond tuples (element 1, element 2, bond length)
        """
        ins_file.add_no_arg_command("ACTA")
        self.driver.run_SHELXTL(ins_file)
        ins_file.remove_command("ACTA")

        with open(self.driver.cif_file) as f:
            cif_file = CifParser.from_string(f.read())

        cif_dict = cif_file.as_dict().values()[0]
        return zip(cif_dict["_geom_bond_atom_site_label_1"], cif_dict["_geom_bond_atom_site_label_2"],
                   [float(x.replace("(", "").replace(")", "")) for x in cif_dict["_geom_bond_distance"]])

    def identify_sites(self):
        """
        If the initial r1 is large, use bond lengths instead of the r1 score as the criteria for identifying which
        electron density peaks are atom sites
        """
        ins_file = self.driver.get_res_file()
        shortest_possible_bond = self.get_shortest_bond(ins_file)
        for i in range(5):
            to_delete = set()
            for bond in sorted(self.get_bonds(ins_file), key=lambda tup: tup[2]):
                # Threshold on how short the bonds are
                if bond[2] / shortest_possible_bond < 0.5:
                    a1_num = int(re.search('\d+', bond[0]).group(0))
                    a2_num = int(re.search('\d+', bond[1]).group(0))
                    to_delete.add(max(a1_num, a2_num))

            ins_file.remove_sites_by_number(to_delete)
            self.run_iter(ins_file)
            ins_file = self.driver.get_res_file()

        # Also remove site if displacement is crazy high
        displacements = [x.displacement for x in ins_file.crystal_sites]
        site_numbers = [x.site_number for x in ins_file.crystal_sites]
        for i, d in enumerate(displacements):
            if d > 1.0:
                to_delete.add(site_numbers[i])
        ins_file.remove_sites_by_number(to_delete)


    def switch_elements(self):
        ins_file = self.driver.get_res_file()

        # Want to make changes from largest displacement to smallest
        displacements = zip(map((lambda x: x.displacement), ins_file.crystal_sites), map((lambda x: x.site_number), ins_file.crystal_sites))
        displacements.sort(key=lambda x: x[0], reverse=True)
        order = [x[1] - 1 for x in displacements]

        num_elems = len(ins_file.elements)
        for i in order:
            for elem in range(1, num_elems+1):

                ins_file.change_element(i, elem)
                self.run_iter(ins_file)

                # Check to make sure no negative displacements caused
                res_file = self.driver.get_res_file()
                cs = res_file.crystal_sites[res_file.get_crystal_sites_by_number(i+1)[0]]
                if cs.displacement < 0.0:
                    self.r1_history[-1] = np.infty


            best_elem = np.argmin(self.r1_history[-num_elems:]) + 1
            ins_file.change_element(i, best_elem)
        self.run_iter(ins_file)

    def try_anisotropy(self):
        """
        Test if adding anisotropy reduces R1 value.  If it does, do so.
        :return:
        """

        ins_file = self.driver.get_res_file()
        prev_ins = copy.deepcopy(ins_file)

        #  Try with anisotropy
        ins_file.add_anisotropy()
        self.run_iter(ins_file)

        #  If anisotropy did not help, revert the ins file
        if self.r1_history[-2] < self.r1_history[-1]:
            self.run_iter(prev_ins)

    def try_exti(self):
        """
        Test if adding extinguishing reduces R1 value.  If it does, do so.
        :return:
        """

        ins_file = self.driver.get_res_file()
        prev_ins = copy.deepcopy(ins_file)

        #  Try with extinguishing
        ins_file.add_exti()
        # print self.r1_history

        self.run_iter(ins_file)


        #  If exti did not help, revert the ins file
        if self.r1_history[-2] < self.r1_history[-1]:
            self.run_iter(prev_ins)
        # print self.r1_history

    def try_add_q(self):
        """
        Try adding q peaks to main crystal sites if it decreases R value
        :return:
        """

        r_before = self.r1_history[-1]
        ins_file = self.driver.get_res_file()
        prev_ins = copy.deepcopy(ins_file)

        # This threshold could be scaled based on the potential atoms
        if ins_file.q_peaks[0].electron_density > 10:
            ins_file.move_q_to_crystal()
            # Find best element for new site
            num_elems = len(ins_file.elements)
            for elem in range(1, num_elems + 1):
                ins_file.change_element(len(ins_file.crystal_sites)-1, elem)
                self.run_iter(ins_file)
            best_elem = np.argmin(self.r1_history[-num_elems:]) + 1
            ins_file.change_element(len(ins_file.crystal_sites)-1, best_elem)
            self.run_iter(ins_file)

            #   If adding one peak helped, recursively try adding another peak until it stops helping
            try_another = True
            if not self.driver.check_res_file():
                try_another = False
            else:
                res_file = self.driver.get_res_file()
                displacements = [x.displacement for x in res_file.crystal_sites]
                if displacements[-1] > np.mean(displacements[:-1]) + 2.0*np.std(displacements[:-1]):
                    try_another = False
            if self.r1_history[-1] > r_before:
                try_another = False
            if try_another:
                self.try_add_q()

            #  If adding peak didn't help, take it back off
            else:
                self.run_iter(prev_ins)


    def use_suggested_weights(self):
        """
        Stop re-initializing weights each time--use previously suggested weights
        :return:
        """
        ins_file = self.driver.get_res_file()
        ins_file.remove_command("WGHT")
        ins_file.commands.append(("WGHT", ins_file.suggested_weight_vals))
        self.run_iter(ins_file)

    def try_remove_site(self, use_ml_model=False):
        """
        Remove crystal sites if they result in bond distances that are too short
        :param ml_model: if True, a machine learning model is used to predict the correct bond length.
                         if False, the ideal bond length is approximated as the sum of the atomic radii
        :return:
        """

        ins_file = self.driver.get_res_file()
        prev_ins = copy.deepcopy(ins_file)
        r_before = self.r1_history[-1]
        r_penalty = 1.1
        bonds = self.get_bonds(ins_file)
        threshold = 0.1

        displacements = [x.displacement for x in ins_file.crystal_sites]

        while True:
            ins_file = copy.deepcopy(prev_ins)
            to_delete = set()
            for a1, a2, distance in bonds:
                ideal_distance = utils.get_ideal_bond_length(a1, a2, use_ml_model)
                # if the distance is too small, remove the lower density atom
                if (ideal_distance - distance) / ideal_distance > threshold:
                    a1_num = int(re.search('\d+', a1).group(0))
                    a2_num = int(re.search('\d+', a2).group(0))
                    to_delete.add(max(a1_num, a2_num))
            ins_file.remove_sites_by_number(to_delete)
            self.run_iter(ins_file)


            if self.r1_history[-1] < r_before * r_penalty or len(to_delete) == 0:
                break
            threshold *= 1.1

        # Also remove site if displacement is crazy high
        to_delete = set()
        ins_file = self.driver.get_res_file()
        displacements = [x.displacement for x in ins_file.crystal_sites]
        site_numbers = [x.site_number for x in ins_file.crystal_sites]
        for i, d in enumerate(displacements):
            if d > 0.6:
                to_delete.add(site_numbers[i])
        ins_file.remove_sites_by_number(to_delete)


    def get_shortest_bond(self, ins_file):
        return sorted([utils.get_ideal_bond_length(el.capitalize(), el.capitalize()) for el in ins_file.elements])[0]

    def change_occupancy(self):
        ins_file = self.driver.get_res_file()

        # Want to make changes from largest displacement to smallest
        displacements = map((lambda x: x.displacement), ins_file.crystal_sites)

        for i, displacement in sorted(enumerate(displacements), key=lambda tup: -tup[1]):
            # don't change occupancy of mixed sites
            if ins_file.crystal_sites[i].site_number in ins_file.mixed_site_numbers:
                continue
            ins_file = self.driver.get_res_file()
            prev_ins = copy.deepcopy(ins_file)
            r_before = self.r1_history[-1]

            # only change occupancy if displacement is >= 2 std deviations away from mean
            mean = np.mean(displacements[:i] + displacements[i+1:])
            std = np.std(displacements[:i] + displacements[i+1:])
            if abs(displacement - mean) / std < 2.0:
                break
            ins_file = self.driver.get_res_file()
            ins_file.add_variable_occupancy(i)
            res = self.run_iter(ins_file)
            # if r1 or displacement go up, undo
            if self.r1_history[-1] > r_before or res.crystal_sites[i].displacement > displacement or float(res.fvar_vals[-1]) < 0.0:
                self.run_iter(prev_ins)

    def try_site_mixing(self):
        """
        Try allowing site mixing, where a single crystal site might have mixed occupancy between two different elements
        Only handles pair-wise mixing.
        :param driver: SHELX driver
        """
        ins_file = self.driver.get_res_file()
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

        tried = []
        # Keep adding site mixing until we've already tried adding site mixing at the top priority sites
        while True:
            bonds = self.get_bonds(ins_file)
            mixing_priority = utils.site_mixing_priority(bonds)
            # In case of ties, find all top tied priorities
            top_priority_score = mixing_priority[0][1]
            top_priority = [priority[0] for priority in mixing_priority if np.abs(priority[1] - top_priority_score) < 0.001]

            # For each of these top priorities, try site mixing
            fail = 0
            for i in top_priority:
                if i in tried:
                    fail += 1
                    continue
                tried.append(i)
                local_ins_history = [self.ins_history[-1]]
                local_r1_history = [self.r1_history[-1]]
                prev_ins = copy.deepcopy(ins_file)
                for pair in pairs:
                    ins_file = copy.deepcopy(prev_ins)
                    ins_file.add_site_mixing(site_number=i, mixing_element_indices=pair)
                    self.run_iter(ins_file, ins_history=local_ins_history, r1_history=local_r1_history)
                    occupancy_var = float(self.driver.get_res_file().fvar_vals[-1])

                    # If the occupancy isn't really split, undo it
                    if occupancy_var < 0.02 or occupancy_var > 0.98:
                        local_ins_history = local_ins_history[:-1]
                        local_r1_history = local_r1_history[:-1]

                best_idx = np.argmin(local_r1_history)
                ins_file = local_ins_history[best_idx]
                self.run_iter(ins_file)
            if fail == len(top_priority):
                break


def test_all(path_to_SXTL_dir, ins_folder, input_prefix="absfac1", output_prefix="temp", use_wine=False):
    subdirs = os.listdir(ins_folder)
    for dirname in subdirs:
        if dirname[0] != "." and dirname[0] != "!":
            print dirname
            test_single(path_to_SXTL_dir, os.path.join(ins_folder, dirname), input_prefix, output_prefix, use_wine=use_wine)


def test_single(path_to_SXTL_dir, dirname, input_prefix="absfac1", output_prefix="temp", print_files=False, use_wine=False):
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
        except Exception:
            print "File structure failure"
            print "~" * 50
            return


    opt = run_single(path_to_SXTL_dir, ins_path, input_prefix, output_prefix, use_wine=use_wine)
    opt_r1 = opt.r1_history[-1]

    r1_tol = 2e-4
    anton_r1 = float(re.search("REM R1 =  (\d\.\d+)", open(final_res).read()).group(1))
    print "Initial r1 = {}".format(opt.r1_history[0])
    print "Optimizer r1 = {}".format(opt_r1)
    print "Reference r1 = {}".format(anton_r1)
    if opt_r1 - anton_r1 > r1_tol:
        print "Not success!"
    if print_files:
        print "Optimizer final result:"
        print open(os.path.join(ins_path, output_prefix + ".res")).read()
        print "Reference final result:"
        print open(final_res).read()
    print "~" * 50


def run_single(path_to_SXTL_dir, ins_path, input_prefix="absfac1", output_prefix="temp", use_wine=False):
    opt = Optimizer()
    opt.run(os.path.join(path_to_SXTL_dir, "xl"), os.path.join(path_to_SXTL_dir, "xs"), ins_path, input_prefix, output_prefix, use_wine=use_wine)
    return opt


def run_all(path_to_SXTL_dir, ins_folder, input_prefix="absfac1", output_prefix="temp", use_wine=False):
    subdirs = os.listdir(ins_folder)
    print subdirs
    for dirname in subdirs:
        if dirname[0] != ".":
            print dirname
            opt = run_single(path_to_SXTL_dir, os.path.join(ins_folder, dirname), input_prefix, output_prefix, use_wine)
            print "Initial r1: {}".format(opt.r1_history[0])
            print "Final r1: {}".format(opt.r1_history[-1])

def main():
    # path_to_SXTL_dir = "/Users/eantono/Documents/program_files/xtal_refinement/SXTL/"
    # ins_folder = "/Users/eantono/Documents/project_files/xtal_refinement/copy/"
    # ins_folder = "/Users/eantono/Documents/project_files/xtal_refinement/!UNSEEN 4-2-1-4/"
    # subdir = "Er4Ru2InGe4"
    path_to_SXTL_dir = "/Users/julialing/Documents/GitHub/crystal_refinement/shelxtl/SXTL/"
    ins_folder = "/Users/julialing/Documents/DataScience/crystal_refinement/single_crystal_data/"

    test_all(path_to_SXTL_dir, ins_folder, use_wine=True)
    # run_all(path_to_SXTL_dir, ins_folder, use_wine=True)
    # test_single(path_to_SXTL_dir, os.path.join(ins_folder, subdir + "/"), "absfac1", print_files=True)
    # run_single(path_to_SXTL_dir, os.path.join(ins_path, "work/"), "absfac1")

if __name__ == "__main__":
    main()


