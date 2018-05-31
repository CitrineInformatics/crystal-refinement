import copy

from crystal_refinement.Optimizer import Optimizer
from crystal_refinement.SHELX.SHELXDriver import SHELXDriver

# input_prefix="input"
# path_to_ins="/Users/eantono/Documents/src/xtal_refinement/examples"
input_prefix="temp"
# path_to_ins="/Users/eantono/Documents/project_files/xtal_refinement/testing_021618_A/mixing"
path_to_ins="/Users/eantono/Documents/project_files/xtal_refinement/testing_021618_A/bond length ml"
# output_prefix="output"
path_to_xl="/Users/eantono/Documents/program_files/xtal_refinement/SXTL/xl.exe"
path_to_xs="/Users/eantono/Documents/program_files/xtal_refinement/SXTL/xs.exe"


driver = SHELXDriver(ins_path=path_to_ins, prefix=input_prefix, path_to_xl=path_to_xl, path_to_xs=path_to_xs,
                     use_wine=True, suppress_ouput=True)


ins_file = copy.deepcopy(driver.get_ins_file())
ins_file.add_command("ACTA")
ins_file.remove_command("L.S.")
ins_file.add_command("L.S.", ["1"])
res = driver.run_SHELXTL(ins_file)
quit()

# copy 1.ins to test.ins
# driver.run_SHELXTL_command("xs")
# quit()

res = driver.get_res_file()
for i in reversed(range(1, len(res.crystal_sites))):
    res.move_crystal_to_q(i)

import numpy as np
res.crystal_sites[0].position = np.asarray([0.0, 0.0, 0.0])
res.crystal_sites[0].occupancy_prefix = 1
res.crystal_sites[0].occupancy = 0.02083
res.crystal_sites[0].displacement = 0.05
driver.run_SHELXTL(res, "xl")

quit()
print(res.get_ins_text())

driver.run_SHELXTL(res, 'xl')

res = driver.get_res_file()
res.move_q_to_crystal()
res.change_element(len(res.crystal_sites), 1)

print(res.get_ins_text())
quit()
driver.run_SHELXTL(res, 'xl')

quit()

optimizer = Optimizer(
    input_prefix=input_prefix,
    path_to_ins=path_to_ins,
    output_prefix=output_prefix,
    path_to_xl=path_to_xl,
    path_to_xs=path_to_xs,
    use_wine=True)
#
#
# optimizer.run()

#
# bonds = utils.get_bonds(driver, ins_file=insfile)
# for bond in bonds:
#     print(bond)
# print(utils.score_compound_bonds(bonds, insfile))