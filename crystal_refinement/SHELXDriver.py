import subprocess, tempfile, shutil, os
from crystal_refinement.SHELXFile import SHELXFile


class SHELXDriver:
    """
    Interface that takes SHELXFile objects and does SHELX things (runs xl, xs, reads in and writes out files)
    """
    def __init__(self, ins_path, prefix, path_to_xl, path_to_xs, use_wine=False):
        """
        :param ins_path: Path to initial ins file (output of xprep)
        :param prefix: Prefix of initial ins file
        :param path_to_SXTL_dir: Path to the directory where xl and xs executables are
        :param use_wine: Set this flag to true if running windows executables using wine on a mac.
        """
        self.use_wine = use_wine
        self.path_to_xl = path_to_xl
        self.path_to_xs = path_to_xs
        self.directory = ins_path
        self.file_prefix = prefix
        self.hkl_file = os.path.join(self.directory, prefix + ".hkl")
        self.ins_file = os.path.join(self.directory, prefix + ".ins")
        self.res_file = os.path.join(self.directory, prefix + ".res")
        self.cif_file = os.path.join(self.directory, prefix + ".cif")
        os.chdir(self.directory)

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
        output = self.run_SHELXTL_command(suppress_output=suppress_output, cmd=cmd)
        if "** Absolute structure probably wrong - invert and repeat refinement **" in output:
            with open(self.ins_file, 'w') as f:
                ins_file_obj.add_command("MOVE", values=["1", "1", "1", "-1"])
                f.write(ins_file_obj.get_ins_text())
            self.run_SHELXTL_command(suppress_output=suppress_output, cmd=cmd)
        if not self.has_valid_res_file():
            return None
        return self.get_res_file()

    def run_SHELXTL_command(self, cmd="xl", suppress_output=True):
        """
        Runs the shelx command
        :param cmd: command to call
        :param suppress_output: Do you want it to print out lots of intermediate results?
        """
        if cmd == "xs":
            command_args = [self.path_to_xs, self.file_prefix]
        else:
            command_args = [self.path_to_xl, self.file_prefix]
        if self.use_wine:
            command_args = ["wine"] + command_args
        output = subprocess.check_output(command_args)
        if not suppress_output:
            print(output)
        return output.decode('utf-8')

    def has_valid_res_file(self):
        """
        Check to make sure that results file exists, has nonzero size, and has q peaks (so it ran properly)
        :return:
        """
        if not os.path.isfile(self.res_file):
            return False
        if os.path.getsize(self.res_file) < 1.0:
            return False
        try:
            res_file = self.get_res_file()
        except Exception:
            return False
        if len(res_file.q_peaks) == 0:
            return False
        return True

