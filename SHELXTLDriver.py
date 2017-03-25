import subprocess, tempfile, shutil, os

# flag for operating system
# interface that takes and returns SHELXTLFile objects and does SHELXTL things

class SHELXTLDriver():
    def __init__(self, original_file_prefix, path_to_SXTL_dir, is_macOS=True):
        self.is_macOS = is_macOS
        self.path_to_SXTL_dir = path_to_SXTL_dir
        self.directory = tempfile.mkdtemp()
        self.file_prefix = os.path.join(self.directory + "temp")
        shutil.copy(original_file_prefix + ".hkl", self.file_prefix + ".hkl")
        shutil.copy(original_file_prefix + ".ins", self.file_prefix + ".ins")
        print "ls in temp dir:"
        print os.listdir(self.directory)
        self.run_SHELXTL_command()
        print "file:"
        with open(self.file_prefix + ".res") as f:
            print f.read()
        pass

    # takes .ins file and returns resulting .res file
    def run_SHELXTL(self, ins_file_obj):
        # return SHELXTLFile()
        pass

    def run_SHELXTL_command(self):
        command_args = [os.path.join(self.path_to_SXTL_dir, "xl.exe"), self.file_prefix]
        if self.is_macOS:
            command_args = ["wine"] + command_args
        subprocess.call(command_args)


def main():
    test_path = "/Users/eantono/Documents/project_files/xtal_refinement/example/4"
    path_to_SXTL_dir = "/Users/eantono/Documents/program_files/xtal_refinement/SXTL"
    SHELXTLDriver(test_path, path_to_SXTL_dir, True)


if __name__ == "__main__":
    main()