from SHELXDriver import SHELXDriver
import os, re, copy, math, itertools, random
import numpy as np
import shutil
from pymatgen.io.cif import CifParser
from pymatgen.core.composition import Element
import utils


class Optimizer:
    def __init__(self):
        self.ins_history = []
        self.r1_history = []

    def run(self, path_to_SXTL_dir, ins_path, input_prefix, output_prefix):
        os.chdir(ins_path)
        shutil.copy(ins_path + input_prefix + ".hkl", ins_path + output_prefix + ".hkl")
        shutil.copy(ins_path + input_prefix + ".ins", ins_path + output_prefix + ".ins")
        self.driver = SHELXDriver(ins_path=ins_path, prefix=output_prefix, path_to_SXTL_dir=path_to_SXTL_dir, use_wine=True)

        # Run first iteration
        self.driver.run_SHELXTL_command(cmd="xs")
        shutil.copy(ins_path + output_prefix + ".res", ins_path + output_prefix + ".ins")

        # Read in and run initial SHELXTL file
        ins_file = self.driver.get_ins_file()
        self.run_iter(ins_file)
        # print "#" * 50
        # Optimization
        if self.r1_history[-1] > 0.1:
            self.identify_sites()
        else:
            self.try_add_q()
            self.try_remove_site()
        # quit()
        self.switch_elements()
        # There is currently an inherent assumption that switch elements will never be called after site mixing, since
        # after site mixing, the indices don't line up anymore
        self.try_site_mixing()
        # print self.driver.get_ins_file().filetxt
        # print "#"*50
        # print self.driver.get_res_file().filetxt
        self.change_occupancy()
        self.try_exti()
        self.try_anisotropy()
        self.use_suggested_weights()
        self.use_suggested_weights()

    def run_iter(self, ins_file, ins_history=None, r1_history=None):
        """
        Run the given ins file through SHELXTL and record the file and resulting r1
        :param self.driver:
        :param ins_file:
        :return:
        """
        if ins_history is None:
            ins_history = self.ins_history

        if r1_history is None:
            r1_history = self.r1_history
        ins_history.append(copy.deepcopy(ins_file))
        res = self.driver.run_SHELXTL(ins_file)
        r1_history.append(res.r1)
        return res

    def identify_sites(self):
        ins_file = self.driver.get_res_file()
        shortest_possible_bond = self.get_shortest_bond(ins_file)
        # print shortest_possible_bond
        for i in range(5):
            to_delete = set()
            # print "worst: {}".format(sorted(self.get_bonds(ins_file), key=lambda tup: tup[2])[0]), sorted(self.get_bonds(ins_file), key=lambda tup: tup[2])[0][2]/shortest_possible_bond
            for bond in sorted(self.get_bonds(ins_file), key=lambda tup: tup[2]):
                # Threshold on how short the bonds are
                if bond[2] / shortest_possible_bond < 0.5:
                    # print bond, bond[2]/shortest_possible_bond
                    a1_num = int(re.search('\d+', bond[0]).group(0))
                    a2_num = int(re.search('\d+', bond[1]).group(0))
                    to_delete.add(max(a1_num, a2_num))
            # print "removing sites {}".format(to_delete)
            ins_file.remove_sites_by_number(to_delete)
            self.run_iter(ins_file)
            ins_file = self.driver.get_res_file()

    def is_converged(self, count, max_iter=100):
        converged = False
        if count >= max_iter or self.r1_history[-1] < 0.02:
            converged = True
        return converged

    def switch_elements(self):
        ins_file = self.driver.get_res_file()

        # Want to make changes from largest displacement to smallest
        displacements = map((lambda x: x.displacement), ins_file.crystal_sites)
        order = (np.argsort(np.asarray(displacements))[::-1]).tolist()
        num_elems = len(ins_file.elements)
        for i in order:
            for elem in range(1, num_elems+1):
                ins_file.change_element(i, elem)
                self.run_iter(ins_file)
            print zip(ins_file.elements, self.r1_history[-num_elems:])
            best_elem = np.argmin(self.r1_history[-num_elems:]) + 1
            ins_file.change_element(i, best_elem)
        self.run_iter(ins_file)
        quit()

    def try_anisotropy(self):
        """
        Test if adding anisotropy reduces R1 value.  If it does, do so.
        :param self.driver: SHELXTL self.driver to run ins file
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
        :param self.driver: SHELXTL self.driver to run ins file
        :return:
        """

        ins_file = self.driver.get_res_file()
        prev_ins = copy.deepcopy(ins_file)

        #  Try with extinguishing
        ins_file.add_exti()
        self.run_iter(ins_file)

        #  If exti did not help, revert the ins file
        if self.r1_history[-2] < self.r1_history[-1]:
            self.run_iter(prev_ins)

    def try_add_q(self):
        """
        Try adding q peaks to main crystal sites if it decreases R value
        :param self.driver:
        :return:
        """

        r_before = self.r1_history[-1]
        ins_file = self.driver.get_res_file()
        prev_ins = copy.deepcopy(ins_file)

        # This threshold could be scaled based on the potential atoms
        if ins_file.q_peaks[0].electron_density > 10:
            ins_file.move_q_to_crystal()
            num_elems = len(ins_file.elements)
            for elem in range(1, num_elems + 1):
                ins_file.change_element(len(ins_file.crystal_sites)-1, elem)
                self.run_iter(ins_file)
            best_elem = np.argmin(self.r1_history[-num_elems:]) + 1
            ins_file.change_element(len(ins_file.crystal_sites)-1, best_elem)
            self.run_iter(ins_file)

            #   If adding one peak helped, recursively try adding another peak until it stops helping
            if self.r1_history[-1] < r_before:
                self.try_add_q()

            #  If adding peak didn't help, take it back off
            else:
                self.run_iter(prev_ins)


    def use_suggested_weights(self):
        ins_file = self.driver.get_res_file()
        ins_file.remove_command("WGHT")
        ins_file.commands.append(("WGHT", ins_file.suggested_weight_vals))
        self.run_iter(ins_file)

    def try_remove_site(self):
        """
        Try adding q peaks to main crystal sites if it decreases R value
        :param self.driver:
        :return:
        """

        ins_file = self.driver.get_res_file()
        prev_ins = copy.deepcopy(ins_file)
        r_before = self.r1_history[-1]
        r_penalty = 1.1
        bonds = self.get_bonds(ins_file)

        threshold = 0.1
        while True:
            ins_file = copy.deepcopy(prev_ins)
            to_delete = set()
            for a1, a2, distance in bonds:
                # calculate approxiate ideal bond length
                # this should really be covalent radii
                ideal_distance = utils.get_ideal_bond_length(a1, a2)
                # if a1 == "Ge5" or a2 == "Ge5":
                #     print a1, a2, (ideal_distance - distance) / ideal_distance, threshold
                # if the distance is too small, remove the lower density atom
                if (ideal_distance - distance) / ideal_distance > threshold:
                    a1_num = int(re.search('\d+', a1).group(0))
                    a2_num = int(re.search('\d+', a2).group(0))
                    to_delete.add(max(a1_num, a2_num))
            ins_file.remove_sites_by_number(to_delete)
            self.run_iter(ins_file)
            # print to_delete
            # print self.r1_history[-1], r_before
            # print "!!!!!"
            if self.r1_history[-1] < r_before * r_penalty or len(to_delete) == 0:
                break
            threshold *= 1.1

    def get_shortest_bond(self, ins_file):
        return sorted([utils.get_ideal_bond_length(el.capitalize(), el.capitalize()) for el in ins_file.elements])[0]

    def change_occupancy(self):
        ins_file = self.driver.get_res_file()

        # Want to make changes from largest displacement to smallest
        displacements = map((lambda x: x.displacement), ins_file.crystal_sites)

        for i, displacement in sorted(enumerate(displacements), key=lambda tup: -tup[1]):
            # don't change occupancy of mixed sites cuz you don't know how
            if ins_file.crystal_sites[i].site_number in ins_file.mixed_sites:
                continue
            ins_file = self.driver.get_res_file()
            prev_ins = copy.deepcopy(ins_file)
            r_before = self.r1_history[-1]
            # arbitrary threshold for which occupancies to change
            mean = np.mean(displacements[:i] + displacements[i+1:])
            std = np.std(displacements[:i] + displacements[i+1:])
            if abs(displacement - mean) / std < 2:
                break
            ins_file = self.driver.get_res_file()
            ins_file.add_variable_occupancy(i)
            res = self.run_iter(ins_file)
            # if r1 or displacement go up, undo
            if self.r1_history[-1] > r_before or res.crystal_sites[i].displacement > displacement:
                self.run_iter(prev_ins)

    def try_site_mixing(self):
        ins_file = self.driver.get_res_file()

        element_list = [Element(el.capitalize()) for el in ins_file.elements]
        pairs = []
        probability_threshold = 2E-4
        for i1, i2 in itertools.combinations(range(len(element_list)), 2):
            e1 = element_list[i1]
            e2 = element_list[i2]
            sp = utils.get_substitution_probability(e1, e2)
            if sp > probability_threshold:
                pairs.append(([i1, i2], sp))

        pairs = [tup[0] for tup in sorted(pairs, key=lambda tup: -tup[1])]

        # only doing pairwise mixing for now, use the BronKerbosch algo to find larger fully connected subgraphs
        # do better ordering (use bond lengths for prioritization)
        # for i in range(1, ins_file.crystal_sites[-1].site_number + 1):
        tried = []
        for i in range(5):
            bonds = self.get_bonds(ins_file)
            mixing_priority = utils.site_mixing_priority(bonds)
            i = mixing_priority[0]
            if i in tried:
                i = mixing_priority[1]
                if i in tried:
                    break
            tried.append(i)
            local_ins_history = [self.ins_history[-1]]
            local_r1_history = [self.r1_history[-1]]
            prev_ins = copy.deepcopy(ins_file)
            for pair in pairs:
                ins_file = copy.deepcopy(prev_ins)
                ins_file.add_site_mixing(i, pair)
                self.run_iter(ins_file, ins_history=local_ins_history, r1_history=local_r1_history)
                occupancy_var = float(self.driver.get_res_file().fvar_vals[-1])
                if occupancy_var < 0.0 or occupancy_var > 1.0:
                    local_ins_history = local_ins_history[:-1]
                    local_r1_history = local_r1_history[:-1]

            best_idx = np.argmin(local_r1_history)
            ins_file = local_ins_history[best_idx]
            self.run_iter(ins_file)


    def get_bonds(self, ins_file):
        ins_file.add_no_arg_command("ACTA")
        self.driver.run_SHELXTL(ins_file)
        ins_file.remove_command("ACTA")

        with open(self.driver.cif_file) as f:
            cif_file = CifParser.from_string(f.read())

        cif_dict = cif_file.as_dict().values()[0]
        return zip(cif_dict["_geom_bond_atom_site_label_1"], cif_dict["_geom_bond_atom_site_label_2"],
                [float(x.replace("(", "").replace(")", "")) for x in cif_dict["_geom_bond_distance"]])

def test_all(path_to_SXTL_dir, ins_folder, input_prefix="absfac1", output_prefix="temp"):
    subdirs = os.listdir(ins_folder)
    for dirname in subdirs:
        if dirname[0] != "." and dirname[0] != "!":
            print dirname
            test_single(path_to_SXTL_dir, os.path.join(ins_folder, dirname), input_prefix, output_prefix)


def test_single(path_to_SXTL_dir, dirname, input_prefix="absfac1", output_prefix="temp", print_files=False):
    try:
        ins_path = os.path.join(dirname, "work") + "/"
        for filename in os.listdir(os.path.join(dirname, "Anton")):
            if ".hkl" in filename:
                shutil.copy(os.path.join(dirname, "Anton", filename),
                            os.path.join(ins_path, input_prefix + ".hkl"))
            if ".res" in filename:
                final_res = os.path.join(dirname, "Anton", filename)
    except Exception:
        print "File structure failure"
        print "~" * 50
        return

    try:
        opt = run_single(path_to_SXTL_dir, ins_path, input_prefix, output_prefix)
        opt_r1 = opt.r1_history[-1]
    except Exception, e:
        print "Optimizer failure: {}".format(e)
        print "~" * 50
        return
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


def run_single(path_to_SXTL_dir, ins_path, input_prefix="absfac1", output_prefix="temp"):
    opt = Optimizer()
    opt.run(path_to_SXTL_dir, ins_path, input_prefix, output_prefix)
    return opt


def main():
    path_to_SXTL_dir = "/Users/eantono/Documents/program_files/xtal_refinement/SXTL/"
    ins_folder = "/Users/eantono/Documents/project_files/xtal_refinement/copy/"
    subdir = "Er4Ru2InGe4"
    # path_to_SXTL_dir = "/Users/julialing/Documents/GitHub/crystal_refinement/shelxtl/SXTL/"
    # ins_path = "/Users/julialing/Documents/DataScience/crystal_refinement/temp/"

    # test_all(path_to_SXTL_dir, ins_folder)
    test_single(path_to_SXTL_dir, os.path.join(ins_folder, subdir + "/"), "absfac1", print_files=True)
    # run_single(path_to_SXTL_dir, os.path.join(ins_path, "work/"), "absfac1")

if __name__ == "__main__":
    main()


