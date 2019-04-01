import subprocess

def main():
    print(subprocess.check_output(["crystal_refinement",
                                   "-d", "/Users/eantono/Documents/project_files [C_CITRINE]/xtal_refinement/paper examples/mixing",
                                   "-c", "/Users/eantono/Documents/src/other/xtal_refinement/working.yml",
                                   "-i", "1",
                                   "-o", "temp",
                                   "-r"
                                   ]))



if __name__ == "__main__":
    main()