from crystal_refinement.SHELXFile import SHELXFile
from crystal_refinement.SHELXDriver import SHELXDriver
from crystal_refinement.Optimizer import Optimizer

optimizer = Optimizer(
    input_prefix="input",
    path_to_ins="/Users/eantono/Documents/src/xtal_refinement/examples",
    output_prefix="output",
    path_to_xl="/Users/eantono/Documents/program_files/xtal_refinement/SXTL/xl.exe",
    path_to_xs="/Users/eantono/Documents/program_files/xtal_refinement/SXTL/xs.exe",
    use_wine=True)


optimizer.run()

#
# bonds = utils.get_bonds(driver, ins_file=insfile)
# for bond in bonds:
#     print(bond)
# print(utils.score_compound_bonds(bonds, insfile))