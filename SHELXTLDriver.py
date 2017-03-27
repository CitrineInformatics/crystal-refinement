import subprocess, tempfile, shutil, os
from SHELXTLFile import SHELXTLFile

# flag for operating system
# interface that takes and returns SHELXTLFile objects and does SHELXTL things

class SHELXTLDriver():
    def __init__(self, ins_path, prefix, path_to_SXTL_dir, is_macOS=True):
        self.is_macOS = is_macOS
        self.path_to_SXTL_dir = path_to_SXTL_dir
        self.directory = ins_path
        #self.directory = tempfile.mkdtemp()
        self.file_prefix = os.path.join(self.directory + "temp")
        self.hkl_file = self.file_prefix + ".hkl"
        self.ins_file = self.file_prefix + ".ins"
        self.res_file = self.file_prefix + ".res"
        shutil.copy(ins_path + prefix + ".hkl", self.hkl_file)
        shutil.copy(ins_path + prefix + ".ins", self.ins_file)
        # we may want some special 1st iteration code here that operates on the first round of iteration
        # after XPREP (removing non-atom Q peaks, removing the MOLE lines etc.)
        # self.run_SHELXTL_command()


    def get_ins_file(self):
        with open(self.ins_file) as f:
            ins_text = f.read()
        return SHELXTLFile(ins_text)

    def get_res_file(self):
        with open(self.res_file) as f:
            res_text = f.read()
        return SHELXTLFile(res_text)

    # takes .ins file and returns resulting .res file
    def run_SHELXTL(self, ins_file_obj):
        with open(self.ins_file, 'w') as f:
            f.write(ins_file_obj.write_ins())
        self.run_SHELXTL_command()
        return self.get_res_file()

    def run_SHELXTL_command(self):
        command_args = [os.path.join(self.path_to_SXTL_dir, "xl.exe"), self.file_prefix]
        if self.is_macOS:
            command_args = ["wine"] + command_args
        subprocess.call(command_args)


def main():
    # path_to_SXTL_dir = "/Users/eantono/Documents/program_files/xtal_refinement/SXTL"
    # test_path = "/Users/eantono/Documents/project_files/xtal_refinement/example/temp"
    path_to_SXTL_dir = "/Users/julialing/Documents/GitHub/crystal-refinement/shelxtl/SXTL/"
    ins_path="/Users/julialing/Documents/DataScience/crystal_refinement/temp/"
    prefix = "orig"
    os.chdir(ins_path)
    driver = SHELXTLDriver(ins_path=ins_path, prefix=prefix, path_to_SXTL_dir=path_to_SXTL_dir, is_macOS=True)
    file_obj = driver.get_ins_file()
    file_obj.add_anisotropy()
    file_obj.change_element(1, 1)
    driver.run_SHELXTL(file_obj).write_ins()


if __name__ == "__main__":
    main()