import os
import re
import shutil
import time
from crystal_refinement.SHELXFile import SHELXFile

from crystal_refinement.Optimizer import Optimizer


def test_all(path_to_SXTL_dir, ins_folder, input_prefix="absfac1", output_prefix="temp", use_wine=False, print_files=False,
             generate_graph=False, truncated_graph=False, graph_path=""):
    subdirs = os.listdir(ins_folder)
    # subdirs = ['Y3Ru2Ge3_lt', 'Y3Ru2Ge3_rt_2']

    for dirname in subdirs:
        if dirname[0] != "." and dirname != "temp":
            print(dirname)
            test_single(path_to_SXTL_dir, os.path.join(ins_folder, dirname), input_prefix, output_prefix, use_wine,
                        print_files, generate_graph, truncated_graph, graph_path)


def test_single(path_to_SXTL_dir, dirname, input_prefix="absfac1", output_prefix="temp", use_wine=False, print_files=False,
                generate_graph=False, truncated_graph=False, graph_path="", result_filename=None):
    if "INS-HKL" in os.listdir(dirname):
        final_res = ""
        ins_path = os.path.join(dirname, "INS-HKL")
        if result_filename is not None:
            final_res = os.path.join(ins_path, result_filename)
            input_prefix = input_prefix
        else:
            for filename in sorted(os.listdir(ins_path), key=lambda name: os.path.getctime(os.path.join(ins_path, name))):
                if ".res" in filename:
                    final_res = os.path.join(ins_path, filename)
                    break
            input_prefix = os.path.basename(final_res).split(".")[0]
        graph_name = os.path.basename(dirname)
        ins_from_result(ins_path, result_file=os.path.basename(final_res), input_prefix=input_prefix)
    elif input_prefix == "detect":
        final_res = ""
        ins_path = dirname
        if "result.res" in os.listdir(ins_path):
            final_res = os.path.join(ins_path, "result.res")
        else:
            for filename in sorted(os.listdir(ins_path), key=lambda name: os.path.getctime(os.path.join(ins_path, name))):
                if ".res" in filename:
                    final_res = os.path.join(ins_path, filename)
                    break
        if "result.res" in final_res:
            input_prefix = "1"
        else:
            input_prefix = os.path.basename(final_res).split(".")[0]
        if "result.res" in os.listdir(ins_path):
            final_res = os.path.join(ins_path, "result.res")
        graph_name = os.path.basename(dirname)
        if input_prefix == "":
            return
        print(ins_path, input_prefix, final_res)
        ins_from_result(ins_path, result_file=os.path.basename(final_res), input_prefix=input_prefix)
    else:
        try:
            ins_path = os.path.join(dirname, "work")
            graph_name = os.path.basename(dirname)
            for filename in os.listdir(os.path.join(dirname, "Anton")):
                if ".hkl" in filename:
                    shutil.copy(os.path.join(dirname, "Anton", filename),
                                os.path.join(ins_path, input_prefix + ".hkl"))
                if ".res" in filename:
                    final_res = os.path.join(dirname, "Anton", filename)
        except Exception:
            try:
                ins_path = dirname
                graph_name = os.path.basename(dirname)
                open(os.path.join(dirname, "1.hkl"))
                open(os.path.join(dirname, "1.ins"))
                final_res = os.path.join(dirname, "result.res")
                open(final_res)
            except Exception, e:
                print("File structure failure", e)
                print("~" * 50)
                return

    start = time.time()

    opt = run_single(path_to_SXTL_dir, ins_path, input_prefix, output_prefix, use_wine,
                     generate_graph, truncated_graph, graph_path, graph_name)
    runtime = time.time() - start

    best_history = opt.history.get_best_history()

    r1_tol = 2e-4
    anton_r1 = float(re.search("REM R1 =  (\d\.\d+)", open(final_res).read()).group(1))
    nominal_formula = "".join(["{}{}".format(el.capitalize(), st) for el, st in zip(best_history[-1].ins_file.elements, best_history[-1].ins_file.formula_units)])
    nominal_elements = set(best_history[-1].ins_file.elements)
    result_formula = best_history[-1].res_file.get_analytic_formula()
    result_elements = set()
    for site in best_history[-1].res_file.crystal_sites:
        result_elements.add(site.el_string.upper())
    anton_res_nsites, anton_res_el_lines = test_result(open(final_res).read())
    opt_res_nsites, opt_res_el_lines = test_result(best_history[-1].res_file.get_ins_text())
    print(len(opt.history.leaves), "path(s) tried")
    print("Initial r1 = {}".format(best_history[0].r1))
    print("Optimizer r1 = {}".format(best_history[-1].r1))
    print("Reference r1 = {}".format(anton_r1))
    print("Nominal formula = {}, result formula = {}".format(nominal_formula, result_formula))
    print("Missing elements = {}".format(nominal_elements.difference(result_elements)))
    print("Reference n_sites = {}, Optimizer n_sites = {}".format(anton_res_nsites, opt_res_nsites))
    if anton_res_nsites != opt_res_nsites:
        print("Mismatched number of sites!!")
    print("Reference sites:")
    for line in anton_res_el_lines:
        print(line)
    print("Optimizer sites:")
    for line in opt_res_el_lines:
        print(line)
    print("Runtime = {}s".format(runtime))
    if best_history[-1].r1 - anton_r1 > r1_tol:
        print("Not success!")
    if print_files:
        print("Optimizer final result:")
        print(open(os.path.join(ins_path, output_prefix + ".res")).read())
        print("Reference final result:")
        print(open(final_res).read())
    print("~" * 50)


def run_single(path_to_SXTL_dir, ins_path, input_prefix="absfac1", output_prefix="temp", use_wine=False,
               generate_graph=False, truncated_graph=False, graph_path="", graph_name="out"):
    opt = Optimizer(os.path.join(path_to_SXTL_dir, "xl.exe"), os.path.join(path_to_SXTL_dir, "xs.exe"), ins_path, input_prefix,
            output_prefix, use_wine=use_wine, n_results=5, use_ml_model=False, max_n_leaves=20,
                    score_weighting=0.5, log_output=False)#, mixing_pairs=[("Ga", "Ge")])
    opt.run()
    if generate_graph:
        if truncated_graph:
            opt.history.head.update_dead_branches()
            opt.history.head.generate_truncated_graph(os.path.join(graph_path, graph_name + "_trunc"))
        else:
            opt.history.head.generate_graph(os.path.join(graph_path, graph_name))
    return opt


def run_all(path_to_SXTL_dir, ins_folder, input_prefix="absfac1", output_prefix="temp", use_wine=False,
            generate_graph=False, truncated_graph=False, graph_path=""):
    subdirs = os.listdir(ins_folder)
    for dirname in subdirs:
        if dirname[0] != ".":
            opt = run_single(path_to_SXTL_dir, os.path.join(ins_folder, dirname), input_prefix, output_prefix, use_wine,
                             generate_graph, truncated_graph, graph_path)
            best_history = opt.history.get_best_history()
            print(len(opt.history.leaves), "path(s) tried")
            print("Initial r1: {}".format(best_history[0].r1))
            print("Final r1: {}".format(best_history[-1].r1))


def test_result(res_text):
    m = re.search("SFAC[^\n]+", res_text)
    elements = m.group(0).split()[1:]
    print(elements)
    el_lines = []
    site_coordinates = set()
    for line in res_text.split("\n"):
        split = line.split()
        if len(split) > 1:
            for el in elements:
                if re.match("^" + el + "\d*$", split[0]) is not None:
                    el_lines.append(split[0:5])
                    site_coordinates.add(" ".join(split[2:4]))

    return len(site_coordinates), sorted(el_lines, key=lambda tup: tup[0])

def ins_from_result(folder_path, result_file="result.res", input_prefix="1"):
    from collections import defaultdict
    lines = defaultdict(lambda: [])
    lines["TREF"].append("\n")


    with open(os.path.join(folder_path, result_file)) as f:
        for line in f:
            split = line.split()
            if len(split) == 1:
                lines[split[0]].append("\n")
            elif len(split) > 1:
                lines[split[0]].append(line[len(split[0]):])


    keys = ["TITL", "CELL", "ZERR", "LATT", "SYMM", "SFAC", "UNIT", "TEMP", "SIZE", "TREF", "HKLF", "END"]
    with open(os.path.join(folder_path, input_prefix + ".ins"), 'w') as f:
        for key in keys:
            values = lines[key]
            for val in values:
                f.write("{} {}".format(key, val))



def main():
    path_to_SXTL_dir = "/Users/eantono/Documents/program_files/xtal_refinement/SXTL/"
    ins_folder = "/Users/eantono/Documents/project_files/xtal_refinement/4-2-1-4 INS and HKL files"
    # ins_folder = "/Users/eantono/Documents/project_files/xtal_refinement/!UNSEEN 4-2-1-4/"
    # ins_folder = "/Users/eantono/Documents/project_files/xtal_refinement/Organized_data1/EASY"
    # ins_folder = "/Users/eantono/Documents/project_files/xtal_refinement/Organized_data1/MIXING"
    # ins_folder = "/Users/eantono/Documents/project_files/xtal_refinement/Organized_data2/EASY"
    # ins_folder = "/Users/eantono/Documents/project_files/xtal_refinement/Organized_data2/Partial occupancy"
    # ins_folder = "/Users/eantono/Documents/project_files/xtal_refinement/Organized_data2/mixing"
    # ins_folder = "/Users/eantono/Documents/project_files/xtal_refinement/testing_020418"
    # ins_folder = "/Users/eantono/Documents/project_files/xtal_refinement/testing_020418 bond score"
    subdir = "Nd4Mn2CdGe4crystal"
    graph_output_path = "/Users/eantono/Documents/src/xtal_refinement/output"
    # path_to_SXTL_dir = "/Users/julialing/Documents/GitHub/crystal_refinement/shelxtl/SXTL/"
    # ins_folder = "/Users/julialing/Documents/DataScience/crystal_refinement/single_crystal_data/"
    # from shutil import copyfile
    # for subdir in os.listdir(ins_folder):
    #     if subdir[0] != ".":
    #         folder = os.path.join(ins_folder, subdir, "INS-HKL")
    #         copyfile(os.path.join(folder, "1.hkl"), os.path.join(folder, "result.hkl"))
    # quit()
    test_all(path_to_SXTL_dir, ins_folder, input_prefix="1", use_wine=True, print_files=False,
      generate_graph=True, truncated_graph=True, graph_path=graph_output_path)
    # test_single(path_to_SXTL_dir, os.path.join(ins_folder, subdir), "detect", use_wine=True, print_files=False,
    #   generate_graph=True, truncated_graph=False, graph_path=graph_output_path, result_filename=None)

    # run_all(path_to_SXTL_dir, ins_folder, input_prefix="1", use_wine=True,
    #   generate_graph=True, graph_path=graph_output_path)
    # run_single(path_to_SXTL_dir, os.path.join(ins_folder, subdir), "1", use_wine=True,
    #   generate_graph=True, graph_path=graph_output_path)

if __name__ == "__main__":
    main()