from SHELXTLFile import SHELXTLFile
from SHELXTLDriver import SHELXTLDriver
import os
import numpy as np
import shutil

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
        self.ins_history.append(ins_file)
        res = driver.run_SHELXTL(ins_file)
        res.write_ins()
        self.r1_history.append(res.r1)
        self.switch_elements(driver)
        self.try_add_q(driver)
        self.try_exti(driver)
        self.try_anisotropy(driver)

    def is_converged(self, count, max_iter=100):
        converged = False
        if count >= max_iter or self.r1_history[-1] < 0.02:
            converged = True
        return converged

    def switch_elements(self, driver):
        ins_file = driver.get_ins_file()

        # Want to make changes from largest displacement to smallest
        displacements = map((lambda x: x.displacement), ins_file.crystal_sites)
        order = (np.argsort(np.asarray(displacements)[::-1])).tolist()
        num_elems = len(ins_file.elements)
        for i in order:
            for elem in range(1, num_elems+1):
                ins_file.change_element(i, elem)
                self.ins_history.append(ins_file)
                res = driver.run_SHELXTL(ins_file)
                res.write_ins()
                self.r1_history.append(res.r1)
            best_elem = np.argmin(self.r1_history[-num_elems:]) + 1
            ins_file.change_element(i, best_elem)
        self.ins_history.append(ins_file)
        res = driver.run_SHELXTL(ins_file)
        res.write_ins()
        self.r1_history.append(res.r1)




    def try_anisotropy(self, driver):
        """
        Test if adding anisotropy reduces R1 value.  If it does, do so.
        :param driver: SHELXTL driver to run ins file
        :return:
        """

        #  Try with anisotropy
        ins_file = driver.get_ins_file()
        ins_file.add_anisotropy()
        self.ins_history.append(ins_file)
        res = driver.run_SHELXTL(ins_file)
        res.write_ins()
        self.r1_history.append(res.r1)

        #  Try without anisotropy
        ins_file.remove_anisotropy()
        self.ins_history.append(ins_file)
        res = driver.run_SHELXTL(ins_file)
        res.write_ins()
        self.r1_history.append(res.r1)

        #  If anisotropy helped, re-apply it
        if self.r1_history[-2] < self.r1_history[-1]:
            ins_file.add_anisotropy()
            self.ins_history.append(ins_file)
            res = driver.run_SHELXTL(ins_file)
            res.write_ins()
            self.r1_history.append(res.r1)



    def try_exti(self, driver):
        """
        Test if adding extinguishing reduces R1 value.  If it does, do so.
        :param driver: SHELXTL driver to run ins file
        :return:
        """

        #  Try with extinguishing
        ins_file = driver.get_ins_file()
        ins_file.add_exti()
        self.ins_history.append(ins_file)
        res = driver.run_SHELXTL(ins_file)
        res.write_ins()
        self.r1_history.append(res.r1)

        #  Try without exti
        ins_file.remove_exti()
        self.ins_history.append(ins_file)
        res = driver.run_SHELXTL(ins_file)
        res.write_ins()
        self.r1_history.append(res.r1)

        #  If exti helped, re-apply it
        if self.r1_history[-2] < self.r1_history[-1]:
            ins_file.add_exti()
            self.ins_history.append(ins_file)
            res = driver.run_SHELXTL(ins_file)
            res.write_ins()
            self.r1_history.append(res.r1)

    def try_add_q(self, driver):
        """
        Try adding q peaks to main crystal sites if it decreases R value
        :param driver:
        :return:
        """
        r_before = self.r1_history[-1]
        ins_file = driver.get_res_file()
        ins_file.move_q_to_crystal()
        num_elems = len(ins_file.elements)
        for elem in range(1, num_elems + 1):
            ins_file.change_element(len(ins_file.crystal_sites)-1, elem)
            self.ins_history.append(ins_file)
            res = driver.run_SHELXTL(ins_file)
            res.write_ins()
            self.r1_history.append(res.r1)
        best_elem = np.argmin(self.r1_history[-num_elems:]) + 1
        ins_file.change_element(len(ins_file.crystal_sites)-1, best_elem)
        self.ins_history.append(ins_file)
        res = driver.run_SHELXTL(ins_file)
        res.write_ins()
        self.r1_history.append(res.r1)

        #   If adding one peak helped, recursively try adding another peak until it stops helping
        if res.r1 < r_before:
            self.try_add_q(driver)

        #  If adding peak didn't help, take it back off
        else:
            ins_file.move_crystal_to_q()
            self.ins_history.append(ins_file)
            res = driver.run_SHELXTL(ins_file)
            res.write_ins()
            self.r1_history.append(res.r1)

    def change_occupancy(self):
        pass

    def change_occupancy_prefix(self):
        pass

def main():
    path_to_SXTL_dir = "/Users/julialing/Documents/GitHub/crystal-refinement/shelxtl/SXTL/"
    ins_path = "/Users/julialing/Documents/DataScience/crystal_refinement/temp/"
    input_prefix = "1"
    output_prefix = "temp"
    opt = Optimizer()
    opt.run(path_to_SXTL_dir, ins_path, input_prefix, output_prefix)
    print opt.r1_history


if __name__ == "__main__":
    main()


