import subprocess, tempfile, shutil, os
from SHELXFile import SHELXFile


class SHELXDriver:
    """
    Interface that takes SHELXFile objects and does SHELX things (runs xl, xs, reads in and writes out files)
    """
    def __init__(self, ins_path, prefix, path_to_SXTL_dir, use_wine=False):
        """
        :param ins_path: Path to initial ins file (output of xprep)
        :param prefix: Prefix of initial ins file
        :param path_to_SXTL_dir: Path to the directory where xl and xs executables are
        :param use_wine: Set this flag to true if running windows executables using wine on a mac.
        """
        self.use_wine = use_wine
        self.path_to_SXTL_dir = path_to_SXTL_dir  # Path to the directory where xl and xs executables are
        self.directory = ins_path  # Path to initial ins file
        self.file_prefix = os.path.join(self.directory + "temp")  # This is where the output gets stored
        self.hkl_file = self.directory + prefix + ".hkl"  # These are the inputs
        self.ins_file = self.directory + prefix + ".ins"
        self.res_file = self.directory + prefix + ".res"
        self.cif_file = self.directory + prefix + ".cif"

    def get_ins_file(self):
        """
        Reads in ins file
        :return: SHELXFile object
        """
        with open(self.ins_file) as f:
            ins_text = f.read()
        return SHELXFile(ins_text)

    def get_res_file(self):
        """
        Reads in res file
        :return: SHELXFILE object
        """
        with open(self.res_file) as f:
            res_text = f.read()
        return SHELXFile(res_text)

    def run_SHELXTL(self, ins_file_obj, suppress_output=True, cmd="xl.exe"):
        """
        Takes SHELXFile object, runs cmd, then returns resulting SHELXFile object
        :param ins_file_obj: File object from ins file
        :param suppress_output: Do you want it to print out lots of intermediate results?
        :param cmd: command to run
        :return: File object from res file
        """
        with open(self.ins_file, 'w') as f:
            f.write(ins_file_obj.get_ins_text())
        self.run_SHELXTL_command(suppress_output=suppress_output, cmd=cmd)
        return self.get_res_file()

    def run_SHELXTL_command(self, cmd="xl.exe", suppress_output=False):
        """
        Runs the shelx command
        :param cmd: command to call
        :param suppress_output: Do you want it to print out lots of intermediate results?
        """
        command_args = [os.path.join(self.path_to_SXTL_dir, cmd), self.file_prefix]
        if self.use_wine:
            command_args = ["wine"] + command_args
        if suppress_output:
            subprocess.call(command_args, stdout=open(os.devnull, "w"))
        else:
            subprocess.call(command_args)


def main():
    # path_to_SXTL_dir = "/Users/eantono/Documents/program_files/xtal_refinement/SXTL"
    # test_path = "/Users/eantono/Documents/project_files/xtal_refinement/example/temp"
    path_to_SXTL_dir = "/Users/julialing/Documents/GitHub/crystal_refinement/shelxtl/SXTL/"
    ins_path="/Users/julialing/Documents/DataScience/crystal_refinement/temp/"
    prefix = "orig"
    os.chdir(ins_path)
    driver = SHELXDriver(ins_path=ins_path, prefix=prefix, path_to_SXTL_dir=path_to_SXTL_dir, use_wine=True)
    file_obj = driver.get_ins_file()
    file_obj.add_anisotropy()
    file_obj.change_element(1, 1)
    driver.run_SHELXTL(file_obj).get_ins_text()


if __name__ == "__main__":
    main()