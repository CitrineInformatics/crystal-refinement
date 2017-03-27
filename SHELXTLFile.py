import re
from CrystalSite import CrystalSite

# Parser for .ins and .res SHELLXTL files
class SHELXTLFile():
    def __init__(self, filetxt):
        # Text that will not be modified, which just needs to be stored in order to re-write the file
        self.extra_text = ["", "", ""]
        self.extra_text_section = 0

        self.elements = []
        self.formula_units = []

        # Stores commands as dict.  If command has no value, value is set to None
        self.commands = {}

        self.fvar_vals = []

        # object to store crystal site information
        # ex. SM1   4    0.086929    0.000000    0.881287    10.50000    0.00621    0.00453 =
        #           0.00623    0.00000    0.00199    0.00000
        self.crystal_sites = []

        # crystal site info (for easier distance calculation and translation between atom sites and q sites)
        self.q_peaks = []

        self.r1 = 0.0
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
        starting_element_keys = ["{}1".format(el) for el in self.elements]
        while True:
            line = lines[line_idx]

            if re.match("^\s*$", line) is None:
                split = line.split()
                key = split[0]
                if len(split) == 1:
                    self.commands[key] = None

                else:
                    if key == "FVAR":
                        self.fvar_vals = split[1:]
                    elif key in starting_element_keys:
                        break
                    else:
                        self.commands[key] = split[1:]
            line_idx += 1

        # Crystal site section
        while True:
            line = lines[line_idx]
            if len(line.split()) == 0:
                break
            if "HKLF" in line:
                break
            if "=" in line:
                crystal_info = line.split()[:-1] + lines[line_idx + 1].split()
                line_idx += 1
            else:
                crystal_info = line.split()
            self.crystal_sites.append(CrystalSite(crystal_info))
            line_idx += 1

        # End text section
        self.extra_text_section = 2
        while line_idx < len(lines):
            line = lines[line_idx]
            if "REM R1 =" in line:
                self.r1 = float(re.search("REM R1 =\s*(\d*\.\d+)", line).group(1))
            self.extra_text[self.extra_text_section] += line + "\n"
            if len(line) > 0 and line[0] == "Q" and line[1].isdigit():
                crystal_info = line.split()
                self.q_peaks.append(crystal_info)
            line_idx += 1


    # write ins file (with commands?)
    def write_ins(self):
        res = self.extra_text[0]
        res += "SFAC " + " ".join(self.elements) + "\n"
        res += "UNIT " + " ".join(self.formula_units) + "\n"
        res += self.extra_text[1] + "\n \n"
        for key in self.commands.keys():
            res += key + "  "
            if self.commands[key] is not None:
                res += " ".join(self.commands[key])
            res += "\n"
        res += "FVAR " + " ".join(self.fvar_vals) + "\n"
        res += "\n".join([" ".join(cs.write_line()) for cs in self.crystal_sites]) + "\n"
        res += self.extra_text[2]
        return res


    # various editing methods ...

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