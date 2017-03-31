from SHELXTLFile import SHELXTLFile
from SHELXTLDriver import SHELXTLDriver
import os, re, copy, math
import numpy as np
import shutil
from pymatgen.io.cif import CifParser
from pymatgen.core.composition import Element
from collections import defaultdict


class Optimizer():
    def __init__(self):
        self.ins_history = []
        self.r1_history = []

    def run(self, path_to_SXTL_dir, ins_path, input_prefix, output_prefix):
        os.chdir(ins_path)
        shutil.copy(ins_path + input_prefix + ".hkl", ins_path + output_prefix + ".hkl")
        shutil.copy(ins_path + input_prefix + ".ins", ins_path + output_prefix + ".ins")
        driver = SHELXTLDriver(ins_path=ins_path, prefix=output_prefix, path_to_SXTL_dir=path_to_SXTL_dir, is_macOS=True)

        # Run first iteration
        driver.run_SHELXTL_command(cmd="xs")
        shutil.copy(ins_path + output_prefix + ".res", ins_path + output_prefix + ".ins")

        # Read in and run initial SHELXTL file
        ins_file = driver.get_ins_file()
        self.run_iter(driver, ins_file)

        # Optimization
        self.try_add_q(driver)
        self.try_remove_site(driver)
        self.switch_elements(driver)
        self.change_occupancy(driver)
        # self.try_split_occupancy(driver)
        # quit()
        self.try_exti(driver)
        self.try_anisotropy(driver)

    def run_iter(self, driver, ins_file):
        """
        Run the given ins file through SHELXTL and record the file and resulting r1
        :param driver:
        :param ins_file:
        :return:
        """
        self.ins_history.append(copy.copy(ins_file))
        res = driver.run_SHELXTL(ins_file)
        self.r1_history.append(res.r1)
        return res

    def is_converged(self, count, max_iter=100):
        converged = False
        if count >= max_iter or self.r1_history[-1] < 0.02:
            converged = True
        return converged

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
        :param driver: SHELXTL driver to run ins file
        :return:
        """

        ins_file = driver.get_res_file()
        prev_ins = copy.copy(ins_file)

        #  Try with anisotropy
        ins_file.add_anisotropy()
        self.run_iter(driver, ins_file)

        #  If anisotropy did not help, revert the ins file
        if self.r1_history[-2] < self.r1_history[-1]:
            self.run_iter(driver, prev_ins)

    def try_exti(self, driver):
        """
        Test if adding extinguishing reduces R1 value.  If it does, do so.
        :param driver: SHELXTL driver to run ins file
        :return:
        """

        ins_file = driver.get_res_file()
        prev_ins = copy.copy(ins_file)

        #  Try with extinguishing
        ins_file.add_exti()
        self.run_iter(driver, ins_file)

        #  If exti did not help, revert the ins file
        if self.r1_history[-2] < self.r1_history[-1]:
            self.run_iter(driver, prev_ins)

    def try_add_q(self, driver):
        """
        Try adding q peaks to main crystal sites if it decreases R value
        :param driver:
        :return:
        """

        r_before = self.r1_history[-1]
        ins_file = driver.get_res_file()
        prev_ins = copy.copy(ins_file)

        threshold_distance = 2.0  # No sites allowed within 2 Angstroms of other sites
        if ins_file.q_peaks[0].calc_min_distance_to_others(ins_file.crystal_sites) > threshold_distance:
            ins_file.move_q_to_crystal()
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


    # def try_remove_site_old(self, driver):
    #     """
    #     Try adding q peaks to main crystal sites if it decreases R value
    #     :param driver:
    #     :return:
    #     """
    #     ins_file = driver.get_ins_file()
    #     threshold_distance = 2.0  # No sites allowed within 2 Angstroms of other sites
    #     for i, site in reversed(list(enumerate(ins_file.crystal_sites))):
    #         if site.calc_min_distance_to_others(ins_file.crystal_sites) < threshold_distance:
    #             r_before = self.r1_history[-1]
    #             del ins_file.crystal_sites[i]
    #             self.add_to_ins_history(ins_file)
    #             res = driver.run_SHELXTL(ins_file)
    #             self.r1_history.append(res.r1)
    #
    #             #  If removing peak didn't help, add it back on
    #             if res.r1 > r_before:
    #                 ins_file.crystal_sites.insert(i, site)
    #                 self.add_to_ins_history(ins_file)
    #                 res = driver.run_SHELXTL(ins_file)
    #
    #                 self.r1_history.append(res.r1)
    #             res.write_ins()

    def try_remove_site(self, driver):
        """
        Try adding q peaks to main crystal sites if it decreases R value
        :param driver:
        :return:
        """

        ins_file = driver.get_res_file()
        prev_ins = copy.copy(ins_file)
        r_before = self.r1_history[-1]
        r_penalty = 1.1
        ins_file.commands["ACTA"] = None
        driver.run_SHELXTL(ins_file)

        with open(driver.cif_file) as f:
            cif_file = CifParser.from_string(f.read())

        cif_dict = cif_file.as_dict().values()[0]
        bonds = zip(cif_dict["_geom_bond_atom_site_label_1"], cif_dict["_geom_bond_atom_site_label_2"],
                    cif_dict["_geom_bond_distance"])

        threshold = 0.1
        while True:
            ins_file = copy.copy(prev_ins)
            to_delete = set()
            for a1, a2, distance in bonds:
                distance = float(distance.replace("(", "").replace(")", ""))
                el1 = Element(re.sub('\d', "", a1))
                el2 = Element(re.sub('\d', "", a2))
                # calculate approxiate ideal bond length
                # this should really be covalent radii
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
            ins_file = driver.get_res_file()
            prev_ins = copy.copy(ins_file)
            r_before = self.r1_history[-1]
            # arbitrary threshold for which occupancies to change
            mean = np.mean(displacements[:i] + displacements[i+1:])
            std = np.std(displacements[:i] + displacements[i+1:])
            if abs(displacement - mean) / std < 2:
                break
            ins_file = driver.get_res_file()
            ins_file.add_variable_occupancy(i)
            res = self.run_iter(driver, ins_file)
            # if r1 or displacement go up, undo
            if self.r1_history[-1] > r_before or res.crystal_sites[i].displacement > displacement:
                self.run_iter(driver, prev_ins)


    def try_split_occupancy(self, driver):
        ins_file = driver.get_res_file()
        prev_ins = copy.copy(ins_file)
        r_before = self.r1_history[-1]
        r_penalty = 1.1
        ins_file.commands["ACTA"] = None
        driver.run_SHELXTL(ins_file)

        with open(driver.cif_file) as f:
            cif_file = CifParser.from_string(f.read())

        cif_dict = cif_file.as_dict().values()[0]
        bonds = zip(cif_dict["_geom_bond_atom_site_label_1"], cif_dict["_geom_bond_atom_site_label_2"],
                    cif_dict["_geom_bond_distance"])
        for site in ins_file.crystal_sites:
            coord = 0
            for a1, a2, distance in bonds:
                if a1 == site.name.capitalize() or a2 == site.name.capitalize():
                    coord += 1
            print site.name, coord
        # quit()
        bonds_with_score = []
        for a1, a2, distance in bonds:
            distance = float(distance.replace("(", "").replace(")", ""))
            el1 = Element(re.sub('\d', "", a1))
            el2 = Element(re.sub('\d', "", a2))
            # calculate approxiate ideal bond length
            # this should really be covalent radii
            ideal_distance = (el1.atomic_radius + el2.atomic_radius)
            bonds_with_score.append((a1, a2, distance, ideal_distance, ideal_distance - distance))

        for x in sorted(bonds_with_score, key= lambda tup: -abs(tup[4])):
            print x
        quit()



def main():
    path_to_SXTL_dir = "/Users/eantono/Documents/program_files/xtal_refinement/SXTL/"
    ins_path = "/Users/eantono/Documents/project_files/xtal_refinement/other_example/"
    input_prefix = "1"
    output_prefix = "temp"

    # path_to_SXTL_dir = "/Users/julialing/Documents/GitHub/crystal-refinement/shelxtl/SXTL/"
    # ins_path = "/Users/julialing/Documents/DataScience/crystal_refinement/temp/"
    # input_prefix = "absfac1"
    # output_prefix = "temp"

    opt = Optimizer()
    opt.run(path_to_SXTL_dir, ins_path, input_prefix, output_prefix)
    print opt.r1_history


if __name__ == "__main__":
    main()


