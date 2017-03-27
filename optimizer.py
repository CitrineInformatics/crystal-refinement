from SHELXTLFile import SHELXTLFile
from SHELXTLDriver import SHELXTLDriver
import os

class Optimizer():
    def __init__(self):
        self.ins_history = []
        self.r1_history = []

    def run(self, path_to_SXTL_dir, ins_path, prefix):
        os.chdir(ins_path)
        driver = SHELXTLDriver(ins_path=ins_path, prefix=prefix, path_to_SXTL_dir=path_to_SXTL_dir, is_macOS=True)

        # Read in and run initial SHELXTL file
        ins_file = driver.get_ins_file()
        self.ins_history.append(ins_file)
        res = driver.run_SHELXTL(ins_file)
        res.write_ins()
        self.r1_history.append(res.r1)


    def is_converged(self, count, max_iter=100):
        converged = False
        if count >= max_iter or self.r1_history[-1] < 0.02:
            converged = True
        return converged

    def switch_elements(self, driver):
        ins_file = driver.get_ins_file()
        self.ins_history.append(ins_file)
        res = driver.run_SHELXTL(ins_file)
        res.write_ins()
        self.r1_history.append(res.r1)
        displacements = 
        pass

    @staticmethod
    def try_anisotropy(ins_file):
        pass

    @staticmethod
    def try_exti(ins_file):
        pass

    def add_q(self):
        pass

    def change_occupancy(self):
        pass

def main():
    path_to_SXTL_dir = "/Users/julialing/Documents/GitHub/crystal-refinement/shelxtl/SXTL/"
    ins_path = "/Users/julialing/Documents/DataScience/crystal_refinement/temp/"
    prefix = "orig"
    opt = Optimizer()
    opt.run(path_to_SXTL_dir, ins_path, prefix)
    print opt.r1_history


if __name__ == "__main__":
    main()


