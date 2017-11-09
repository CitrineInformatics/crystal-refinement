from crystal_refinement.SHELXFile import SHELXFile
from crystal_refinement.SHELXDriver import SHELXDriver
import utils

with open("/Users/eantono/Documents/project_files/xtal_refinement/Organized_data2/EASY/mar1229_Rb4Zn7As7/INS-HKL/temp.res") as f:
    insfile = SHELXFile(f.read())

print(insfile.filetxt)

driver = SHELXDriver(ins_path="/Users/eantono/Documents/project_files/xtal_refinement/Organized_data2/EASY/mar1229_Rb4Zn7As7/INS-HKL",
                     prefix="temp",
                     path_to_xl="/Users/eantono/Documents/program_files/xtal_refinement/SXTL/xl.exe",
                     path_to_xs="/Users/eantono/Documents/program_files/xtal_refinement/SXTL/xs.exe",
                      use_wine=True)


bonds = utils.get_bonds(driver, ins_file=insfile)
for bond in bonds:
    print(bond)
print(utils.score_compound_bonds(bonds, insfile))