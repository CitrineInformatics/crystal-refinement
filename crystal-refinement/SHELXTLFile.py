import re, copy
from CrystalSite import CrystalSite

# Parser for .ins and .res SHELLXTL files
class SHELXTLFile():
    def __init__(self, filetxt):
        # Text that will not be modified, which just needs to be stored in order to re-write the file
        self.filetxt = filetxt
        self.extra_text = ["", "", ""]
        self.extra_text_section = 0

        self.elements = []
        self.formula_units = []

        # Stores commands as dict.  If command has no value, value is set to None
        self.commands = []

        self.fvar_vals = []

        # object to store crystal site information
        # ex. SM1   4    0.086929    0.000000    0.881287    10.50000    0.00621    0.00453 =
        #           0.00623    0.00000    0.00199    0.00000
        self.crystal_sites = []
        self.mixed_sites = []

        # crystal site info (for easier distance calculation and translation between atom sites and q sites)
        self.q_peaks = []

        self.r1 = 0.0
        self.suggested_weight_vals = ""
        self.read(filetxt)

    # do actual file parsing
    def read(self, filetxt):
        lines = filetxt.split("\n")
        line_idx = 0

        # First text section, break at SFAC line
        while True:
            line = lines[line_idx]
            if "SFAC" in line:
                break
            self.extra_text[self.extra_text_section] += line + "\n"
            line_idx += 1

        # Parse composition info
        self.elements = lines[line_idx].split()[1:]
        self.formula_units = lines[line_idx + 1].split()[1:]

        line_idx += 2

        # Second text section, break at PLAN line
        self.extra_text_section = 1
        while True:
            line = lines[line_idx]
            self.extra_text[self.extra_text_section] += line + "\n"
            line_idx += 1
            if "PLAN" in line:
                break

        line_idx += 1

        # Command/parameter section break on starting element key
        starting_element_keys = ["{}".format(el) for el in self.elements]
        while True:
            line = lines[line_idx]

            if re.match("^\s*$", line) is None:
                split = line.split()
                key = split[0]
                if key == "FVAR":
                    self.fvar_vals = split[1:]
                elif len(split) == 1:
                    self.commands.append((key, None))
                else:
                    if key == "MOLE":
                        break
                    elif key[:-1] in starting_element_keys:
                        break
                    else:
                        self.commands.append((key, split[1:]))
            line_idx += 1

        # Crystal site section
        while True:
            line = lines[line_idx]
            if len(line.split()) == 0:
                break
            if "HKLF" in line:
                break
            if "MOLE" in line:
                line_idx += 1
                continue
            if "=" in line:
                crystal_info = line.split()[:-1] + lines[line_idx + 1].split()
                line_idx += 1
            else:
                crystal_info = line.split()
            if crystal_info[0][0] == "Q":
                self.q_peaks.append(CrystalSite(crystal_info))
            else:
                self.crystal_sites.append(CrystalSite(crystal_info))
            line_idx += 1

        # End text section
        self.extra_text_section = 2
        while line_idx < len(lines):
            line = lines[line_idx]
            if "REM R1 =" in line:
                self.r1 = float(re.search("REM R1 =\s*(\d*\.\d+)", line).group(1))
            if "WGHT" in line:
                self.suggested_weight_vals = line.split()[1:]
            self.extra_text[self.extra_text_section] += line + "\n"
            if len(line) > 0 and line[0] == "Q" and line[1].isdigit():
                crystal_info = line.split()
                self.q_peaks.append(CrystalSite(crystal_info))
            line_idx += 1


    # write ins file (with commands?)
    def get_ins_text(self):
        res = self.extra_text[0]
        res += "SFAC " + " ".join(self.elements) + "\n"
        res += "UNIT " + " ".join(self.formula_units) + "\n"
        res += self.extra_text[1] + "\n \n"
        for key, values in self.commands:
            res += key + "  "
            if values is not None:
                res += " ".join(values)
            res += "\n"
        res += "FVAR " + " ".join([str(x) for x in self.fvar_vals]) + "\n"
        res += "\n".join([" ".join(cs.write_line()) for cs in self.crystal_sites]) + "\n"
        res += self.extra_text[2]


        return res


    # various editing methods ...

    def add_no_arg_command(self, cmd):
        if cmd not in map(lambda tup: tup[0], self.commands):
            self.commands.append((cmd, None))

    def remove_command(self, cmd):
        cmd_keys = map(lambda tup: tup[0], self.commands)
        if cmd in cmd_keys:
            indices = [i for i, val in enumerate(cmd_keys) if val is cmd]
            indices.reverse()
            for i in indices:
                del self.commands[i]

    def add_anisotropy(self):
        self.add_no_arg_command("ANIS")

    def remove_anisotropy(self):
        self.remove_command("ANIS")

    def add_exti(self):
        self.add_no_arg_command("EXTI")

    def remove_exti(self):
        self.remove_command("EXTI")

    def change_element(self, site_index, element_index):
        self.crystal_sites[site_index].element = element_index
        self.crystal_sites[site_index].name = self.elements[element_index-1] + str(site_index+1)

    def move_q_to_crystal(self):
        self.crystal_sites.append(self.q_peaks[0])
        del self.q_peaks[0]

    def move_crystal_to_q(self, site_idx=-1):
        self.q_peaks.insert(0, self.crystal_sites[site_idx])
        del self.crystal_sites[site_idx]

    def remove_sites_by_number(self, site_numbers):
        self.crystal_sites = [site for site in self.crystal_sites if site.site_number not in site_numbers]

    def add_variable_occupancy(self, site_index):
        site = self.crystal_sites[site_index]
        if site.occupancy_prefix == 1:
            self.fvar_vals.append(0.5)
            site.occupancy_prefix = len(self.fvar_vals)

    def add_site_mixing(self, site_number, mixing_element_indices):
        site_indices = self.get_crystal_sites_by_number(site_number)
        self.commands.append(("EXYZ", [self.elements[i] + str(site_number) for i in mixing_element_indices]))
        self.commands.append(("EADP", [self.elements[i] + str(site_number) for i in mixing_element_indices]))

        # this is 2 site mixing only
        self.fvar_vals.append(0.5)

        mixed_sites = []

        for i, element_idx in enumerate(mixing_element_indices):
            new_site = copy.deepcopy(self.crystal_sites[site_indices[0]])
            new_site.name = self.elements[element_idx] + str(site_number)
            new_site.element = element_idx + 1
            if i == 0:
                new_site.occupancy_prefix = len(self.fvar_vals)
            else:
                new_site.occupancy_prefix = -len(self.fvar_vals)
            mixed_sites.append(new_site)

            # for site in mixed_sites:
            # print site.occupancy_prefix, site.write_line()
        self.remove_sites_by_number([site_number])
        self.crystal_sites.extend(mixed_sites)
        self.crystal_sites.sort(key=lambda site: site.site_number)
        self.mixed_sites.append(site_number)

    def get_crystal_sites_by_number(self, index):
        return [i for i, site in enumerate(self.crystal_sites) if site.site_number == index]


def main():
    test_file = "/Users/julialing/Documents/DataScience/crystal_refinement/4-2-1-4_single_crystal/Example_from_slides/7.ins"

    with open(test_file) as f:
        text = f.read()
        print text
        print "\n\n" + "~"*50 + "\n\n"
        file_obj = SHELXTLFile(text)

        print file_obj.write_ins()


if __name__ == "__main__":
    main()