import os
from crystal_refinement.SHELXDriver import SHELXDriver

def main():
    # path_to_SXTL_dir = "/Users/eantono/Documents/program_files/xtal_refinement/SXTL"
    # test_path = "/Users/eantono/Documents/project_files/xtal_refinement/example/temp"
    path_to_SXTL_dir = "/Users/julialing/Documents/GitHub/crystal_refinement/shelxtl/SXTL/"
    ins_path="/Users/julialing/Documents/DataScience/crystal_refinement/temp/"
    prefix = "orig"
    os.chdir(ins_path)
    driver = SHELXDriver(ins_path=ins_path, prefix=prefix, path_to_xs=path_to_SXTL_dir+"xs", path_to_xl=path_to_SXTL_dir+"xl", use_wine=True)
    file_obj = driver.get_ins_file()
    file_obj.add_anisotropy()
    file_obj.change_element(1, 1)
    driver.run_SHELXTL(file_obj).get_ins_text()

if __name__ == "__main__":
    main()