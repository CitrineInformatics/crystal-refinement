from SHELXFile import SHELXFile
from SHELXDriver import SHELXDriver
import os, re, copy, math, itertools, random
import numpy as np
import shutil
from pymatgen.io.cif import CifParser
from pymatgen.core.composition import Element
from pymatgen.structure_prediction.substitution_probability import SubstitutionProbability
from collections import defaultdict


class Optimizer:
    """
    Class for performing single crystal refinement
    """
    def __init__(self):
        self.ins_history = []  # History of all previous ins files
        self.r1_history = []  # History of all previous R1 values

    def run(self, path_to_SXTL_dir, ins_path, input_prefix, output_prefix, use_wine=False):
        """
        Method to run the optimization
        :param path_to_SXTL_dir: path to xl/xs executables
        :param ins_path: path to ins file output by xprep
        :param input_prefix: prefix of ins file (eg for file.ins, the prefix would be "file")
        :param output_prefix: prefix of the result file that the optimizer will output
        :return:
        """
        # Copy ins and hkl file to output prefix
        os.chdir(ins_path)
        shutil.copy(ins_path + input_prefix + ".hkl", ins_path + output_prefix + ".hkl")
        shutil.copy(ins_path + input_prefix + ".ins", ins_path + output_prefix + ".ins")
        driver = SHELXDriver(ins_path=ins_path, prefix=output_prefix, path_to_SXTL_dir=path_to_SXTL_dir, use_wine=use_wine)

        # Run first iteration using xs
        driver.run_SHELXTL_command(cmd="xs")
        shutil.copy(ins_path + output_prefix + ".res", ins_path + output_prefix + ".ins")

        # Read in and run initial SHELX file
        ins_file = driver.get_ins_file()
        self.run_iter(driver, ins_file)

        # Optimization
        self.try_add_q(driver)
        self.try_remove_site(driver)
        self.switch_elements(driver)
        # There is currently an inherent assumption that switch elements will never be called after site mixing, since
        # after site mixing, the indices don't line up anymore
        self.try_site_mixing(driver)
        self.change_occupancy(driver)
        self.try_exti(driver)
        self.try_anisotropy(driver)
        self.use_suggested_weights(driver)

    def run_iter(self, driver, ins_file, ins_history=None, r1_history=None):
        """
        Run the given ins file through SHELXTL and record the file and resulting r1
        :param driver: SHELXDriver object
        :param ins_file: SHELXFile object
        :return res: SHELXFile object
        """
        if ins_history is None:
            ins_history = self.ins_history

        if r1_history is None:
            r1_history = self.r1_history
        ins_history.append(copy.deepcopy(ins_file))
        res = driver.run_SHELXTL(ins_file)
        r1_history.append(res.r1)
        return res

    def switch_elements(self, driver):
        ins_file = driver.get_res_file()

        # Want to make changes from largest displacement to smallest
        displacements = map((lambda x: x.displacement), ins_file.crystal_sites)
        order = (np.argsort(np.asarray(displacements))[::-1]).tolist()
        num_elems = len(ins_file.elements)
        for i in order:
            for elem in range(1, num_elems+1):
                ins_file.change_element(i, elem)
                self.run_iter(driver, ins_file)
            best_elem = np.argmin(self.r1_history[-num_elems:]) + 1
            ins_file.change_element(i, best_elem)
        self.run_iter(driver, ins_file)

    def try_anisotropy(self, driver):
        """
        Test if adding anisotropy reduces R1 value.  If it does, do so.
        :param driver: SHELX driver to run ins file
        :return:
        """

        ins_file = driver.get_res_file()
        prev_ins = copy.deepcopy(ins_file)

        #  Try with anisotropy
        ins_file.add_anisotropy()
        self.run_iter(driver, ins_file)

        #  If anisotropy did not help, revert the ins file
        if self.r1_history[-2] < self.r1_history[-1]:
            self.run_iter(driver, prev_ins)

    def try_exti(self, driver):
        """
        Test if adding extinguishing reduces R1 value.  If it does, do so.
        :param driver: SHELX driver to run ins file
        :return:
        """

        ins_file = driver.get_res_file()
        prev_ins = copy.deepcopy(ins_file)

        #  Try with extinguishing
        ins_file.add_exti()
        self.run_iter(driver, ins_file)

        #  If exti did not help, revert the ins file
        if self.r1_history[-2] < self.r1_history[-1]:
            self.run_iter(driver, prev_ins)

    def try_add_q(self, driver):
        """
        Try adding q peaks to main crystal sites if it decreases R value
        :param driver: SHELX driver
        :return:
        """

        r_before = self.r1_history[-1]
        ins_file = driver.get_res_file()
        prev_ins = copy.deepcopy(ins_file)
        ins_file.move_q_to_crystal()

        # Find best element for new site
        num_elems = len(ins_file.elements)
        for elem in range(1, num_elems + 1):
            ins_file.change_element(len(ins_file.crystal_sites)-1, elem)
            self.run_iter(driver, ins_file)
        best_elem = np.argmin(self.r1_history[-num_elems:]) + 1
        ins_file.change_element(len(ins_file.crystal_sites)-1, best_elem)
        self.run_iter(driver, ins_file)

        #   If adding one peak helped, recursively try adding another peak until it stops helping
        if self.r1_history[-1] < r_before:
            self.try_add_q(driver)

        #  If adding peak didn't help, take it back off
        else:
            self.run_iter(driver, prev_ins)


    def use_suggested_weights(self, driver):
        """
        Stop re-initializing weights each time--use previously suggested weights
        :param driver:
        :return:
        """
        ins_file = driver.get_res_file()
        ins_file.remove_command("WGHT")
        ins_file.commands.append(("WGHT", ins_file.suggested_weight_vals))
        self.run_iter(driver, ins_file)

    def try_remove_site(self, driver, ml_model=False):
        """
        Remove crystal sites if they result in bond distances that are too short
        :param driver: SHELX driver
        :param ml_model: if True, a machine learning model is used to predict the correct bond length.
                         if False, the ideal bond length is approximated as the sum of the atomic radii
        :return:
        """

        ins_file = driver.get_res_file()
        prev_ins = copy.deepcopy(ins_file)
        r_before = self.r1_history[-1]
        r_penalty = 1.1

        #  Use CIF file to get distances between sites
        ins_file.add_no_arg_command("ACTA")  # To write out CIF file
        driver.run_SHELXTL(ins_file)
        ins_file.remove_command("ACTA")
        with open(driver.cif_file) as f:
            cif_file = CifParser.from_string(f.read())
        cif_dict = cif_file.as_dict().values()[0]
        bonds = zip(cif_dict["_geom_bond_atom_site_label_1"], cif_dict["_geom_bond_atom_site_label_2"],
                    cif_dict["_geom_bond_distance"])

        threshold = 0.1
        while True:
            ins_file = copy.deepcopy(prev_ins)
            to_delete = set()
            for a1, a2, distance in bonds:
                distance = float(distance.replace("(", "").replace(")", ""))  # Take off parenthesis
                el1 = Element(re.sub('\d', "", a1))
                el2 = Element(re.sub('\d', "", a2))
                if not ml_model:
                    # Use pymatgen to get approximate bond length = sum of atomic radii
                    ideal_distance = (el1.atomic_radius + el2.atomic_radius)
                # if the distance is too small, remove the lower density atom
                if (ideal_distance - distance) / ideal_distance > threshold:
                    a1_num = int(re.search('\d+', a1).group(0))
                    a2_num = int(re.search('\d+', a2).group(0))
                    to_delete.add(max(a1_num, a2_num))
            ins_file.remove_sites_by_number(to_delete)
            self.run_iter(driver, ins_file)
            if self.r1_history[-1] < r_before * r_penalty or len(to_delete) == 0:
                break
            threshold *= 1.1

    def change_occupancy(self, driver):
        ins_file = driver.get_res_file()

        # Want to make changes from largest displacement to smallest
        displacements = map((lambda x: x.displacement), ins_file.crystal_sites)

        for i, displacement in sorted(enumerate(displacements), key=lambda tup: -tup[1]):
            # don't change occupancy of mixed sites
            if ins_file.crystal_sites[i].site_number in ins_file.mixed_site_numbers:
                continue
            ins_file = driver.get_res_file()
            prev_ins = copy.deepcopy(ins_file)
            r_before = self.r1_history[-1]

            # only change occupancy if displacement is >= 2 std deviations away from mean
            mean = np.mean(displacements[:i] + displacements[i+1:])
            std = np.std(displacements[:i] + displacements[i+1:])
            if abs(displacement - mean) / std < 2.0:
                break
            ins_file = driver.get_res_file()
            ins_file.add_variable_occupancy(i)
            res = self.run_iter(driver, ins_file)
            # if r1 or displacement go up, undo
            if self.r1_history[-1] > r_before or res.crystal_sites[i].displacement > displacement:
                self.run_iter(driver, prev_ins)

    def get_specie(self, el, ox):
        """
        Return string for ion
        :param el: Element abbreviation (eg. Fe, Au, ...)
        :param ox: Oxidation number (eg. +2)
        :return: eg. Fe2+
        """
        if ox > 0:
            return el.symbol + str(ox) + "+"
        else:
            return el.symbol + str(abs(ox)) + "-"

    def get_substitution_probability(self, el1, el2):
        """
        Use pymatgen to calculate substitution probability for use in site mixing
        :param el1: name of first element
        :param el2: name of second element
        :return: total probability of substitution between two elements
        """
        sp = SubstitutionProbability()
        total = 0

        # Sum over all common oxidation states for both elements
        for o1 in el1.common_oxidation_states:
            for o2 in el2.common_oxidation_states:
                total += sp.prob(self.get_specie(el1, o1), self.get_specie(el2, o2))
        return total

    def try_site_mixing(self, driver):
        """
        Try allowing site mixing, where a single crystal site might have mixed occupancy between two different elements
        Only handles pair-wise mixing.
        :param driver: SHELX driver
        """
        ins_file = driver.get_res_file()

        element_list = [Element(el.capitalize()) for el in ins_file.elements]
        pairs = []
        probability_threshold = 2E-4

        # For all elements in compound, calculate substitution probabilities
        # If substitution probability is > probability_threshold, then save it to pairs list.
        for i1, i2 in itertools.combinations(range(len(element_list)), 2):
            e1 = element_list[i1]
            e2 = element_list[i2]
            sp = self.get_substitution_probability(e1, e2)
            if sp > probability_threshold:
                pairs.append(([i1, i2], sp))

        # Sort pairs by substitution probability (largest to smallest)
        pairs = [tup[0] for tup in sorted(pairs, key=lambda tup: -tup[1])]

        tried = []

        # Keep adding site mixing until we've already tried adding site mixing at the top priority sites
        while True:
            mixing_priority = self.site_mixing_priority(driver, ins_file)
            # In case of ties, find all top tied priorities
            top_priority = [priority for priority in mixing_priority if priority == mixing_priority[0]]

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
                    self.run_iter(driver, ins_file, ins_history=local_ins_history, r1_history=local_r1_history)
                    occupancy_var = float(driver.get_res_file().fvar_vals[-1])

                    # If the occupancy isn't really split, undo it
                    if occupancy_var < 0.02 or occupancy_var > 0.98:
                        local_ins_history = local_ins_history[:-1]
                        local_r1_history = local_r1_history[:-1]

                best_idx = np.argmin(local_r1_history)
                ins_file = local_ins_history[best_idx]
                self.run_iter(driver, ins_file)
            if fail == len(top_priority):
                break


    def site_mixing_priority(self, driver, ins_file):
        """
        Determines priority for adding in site mixing.  Finds most problematic sites based on whether bonds
        are shorter than expected.
        :param driver: SHELX driver
        :param ins_file: SHELX File
        :return:
        """

        # Generate CIF file
        ins_file.add_no_arg_command("ACTA")
        driver.run_SHELXTL(ins_file)
        ins_file.remove_command("ACTA")

        # Use pymatgen to get bond distances
        with open(driver.cif_file) as f:
            cif_file = CifParser.from_string(f.read())

        cif_dict = cif_file.as_dict().values()[0]
        bonds = zip(cif_dict["_geom_bond_atom_site_label_1"], cif_dict["_geom_bond_atom_site_label_2"],
                    cif_dict["_geom_bond_distance"])
        bond_by_atom = defaultdict(lambda: [])
        for a1, a2, distance in bonds:
            distance = float(distance.replace("(", "").replace(")", ""))
            el1 = Element(re.sub('\d', "", a1))
            el2 = Element(re.sub('\d', "", a2))
            # Calculate amount which bonds are shorter than they are supposed to be
            bond_by_atom[a1].append(distance - (el1.atomic_radius + el2.atomic_radius))
            bond_by_atom[a2].append(distance - (el1.atomic_radius + el2.atomic_radius))

        # Sort by which bonds are the shortest compared to what we'd expect
        # Average over deviation from 4 shortest bonds
        sites = sorted(map(lambda tup: (tup[0], sum(sorted(tup[1])[:4])), bond_by_atom.items()), key=lambda tup: tup[1])
        return map(lambda tup: int(re.search("\d+", tup[0]).group(0)), sites)


def main():
    path_to_SXTL_dir = "/Users/eantono/Documents/program_files/xtal_refinement/SXTL/"
    # ins_path = "/Users/eantono/Documents/project_files/xtal_refinement/other_example/"
    ins_folder = "/Users/eantono/Documents/project_files/xtal_refinement/4-2-1-4_single_crystal/"
    # input_prefix = "1"
    output_prefix = "temp"

    # path_to_SXTL_dir = "/Users/julialing/Documents/GitHub/crystal_refinement/shelxtl/SXTL/"
    # ins_path = "/Users/julialing/Documents/DataScience/crystal_refinement/temp/"
    # input_prefix = "absfac1"
    # output_prefix = "temp"
    subdirs = os.listdir(ins_folder)
    # random.shuffle(subdirs)
    for dirname in subdirs:
        if dirname[0] != "." and dirname[0] != "!":
            ins_path = os.path.join(ins_folder, dirname, "work") + "/"
            sorted_files = sorted(os.listdir(ins_path), key=lambda filename: os.path.getmtime(os.path.join(ins_path, filename)))
            input_prefix = ""
            final_res = ""
            for filename in sorted_files:
                if ".hkl" in filename:
                    input_prefix = filename[:-4]
                    if input_prefix + ".ins" in sorted_files:
                        break
            for filename in reversed(sorted_files):
                if output_prefix in filename:
                    continue
                if ".res" in filename:
                    final_res = os.path.join(ins_path, filename)
            print input_prefix
            print final_res
            opt = Optimizer()
            opt.run(path_to_SXTL_dir, ins_path, input_prefix, output_prefix)
            print opt.r1_history
            print open(final_res).read()
            print "#"*50
            print open(os.path.join(ins_path, output_prefix + ".res"))
            quit()

def test_main():
    path_to_SXTL_dir = "/Users/eantono/Documents/program_files/xtal_refinement/SXTL/"
    ins_path = "/Users/eantono/Documents/project_files/xtal_refinement/Ce4Co2InCe4Ex/"
    path_to_SXTL_dir = "/Users/julialing/Documents/GitHub/crystal_refinement/shelxtl/SXTL/"
    ins_path = "/Users/julialing/Documents/DataScience/crystal_refinement/temp/"

    input_prefix = "c2mxs1"
    input_prefix = "absfac1"
    output_prefix = "temp"

    # sorted_files = sorted(os.listdir(ins_path), key=lambda filename: os.path.getmtime(os.path.join(ins_path, filename)))
    # final_res = ""
    # for filename in reversed(sorted_files):
    #     if output_prefix in filename:
    #         continue
    #     if ".res" in filename:
    #         final_res = os.path.join(ins_path, filename)
    #shutil.copy(os.path.join(ins_path, "raw.hkl"), os.path.join(ins_path, input_prefix + ".hkl"))
    opt = Optimizer()
    opt.run(path_to_SXTL_dir, ins_path, input_prefix, output_prefix, use_wine=True)

    print opt.r1_history
    # print open(final_res).read()
    # print "#" * 50
    # print open(os.path.join(ins_path, output_prefix + ".res")).read()


if __name__ == "__main__":
    # main()
    test_main()


