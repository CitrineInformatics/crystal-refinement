from __future__ import absolute_import
import os
import subprocess

from crystal_refinement.SHELX.SHELXFile import SHELXFile


class SHELXDriver:
    """
    Interface that takes SHELXFile objects and does SHELX things (runs xl, xs, reads in and writes out files)
    """
    def __init__(self, directory, prefix, path_to_xl, path_to_xs, use_wine=False, suppress_ouput=True):
        """
        :param directory: Path to initial ins file (output of xprep)
        :param prefix: Prefix of initial ins file
        :param path_to_SXTL_dir: Path to the directory where xl and xs executables are
        :param use_wine: Set this flag to true if running windows executables using wine on a mac.
        """
        self.use_wine = use_wine
        self.path_to_xl = path_to_xl
        self.path_to_xs = path_to_xs
        self.directory = directory
        self.file_prefix = prefix
        self.hkl_file = os.path.join(self.directory, prefix + ".hkl")
        self.ins_file = os.path.join(self.directory, prefix + ".ins")
        self.res_file = os.path.join(self.directory, prefix + ".res")
        self.cif_file = os.path.join(self.directory, prefix + ".cif")
        self.suppress_output = suppress_ouput
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

    def run_SHELXTL(self, ins_file_obj, cmd="xl.exe"):
        """
        Takes SHELXFile object, runs cmd, then returns resulting SHELXFile object
        :param ins_file_obj: File object from ins file
        :param suppress_output: Do you want it to print out lots of intermediate results?
        :param cmd: command to run
        :return: File object from res file
        """
        with open(self.ins_file, 'w') as f:
            f.write(ins_file_obj.to_string())
        output = self.run_SHELXTL_command(cmd=cmd)
        recognized_errors = ["** Cell contents from UNIT instruction and atom list do not agree **",
                             "** Extinction (EXTI) or solvent water (SWAT) correction may be required **"]
        if "** Absolute structure probably wrong - invert and repeat refinement **" in output:
            with open(self.ins_file, 'w') as f:
                ins_file_obj.add_command("MOVE", values=["1", "1", "1", "-1"])
                f.write(ins_file_obj.to_string())
            self.run_SHELXTL_command(cmd=cmd)
        if "** NEGATIVE OCCUPANCY FOR ATOM" in output:
            return None
        if not self.has_valid_res_file():
            return None

        # This currently spams the output too much
        # if "**" in output:
        #     matches = re.findall("\*\*[^(*|\n)]+\*\*", output)
        #     # use a regex, you don't want to throw this out if "any" matches
        #     if not all(x in recognized_errors for x in matches):
        #         print("Unrecognized error in refinement output:")
        #         print(output)
        #         print("INS file:")
        #         print(ins_file_obj.to_string())

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
        if not self.suppress_output or not suppress_output:
            print(output)
        return output

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

