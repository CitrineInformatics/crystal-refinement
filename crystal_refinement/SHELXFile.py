import re, copy
from CrystalSite import CrystalSite
import numpy as np
#TODO JL: REPLACE ALL main() with examples that will run locally


class SHELXFile:
    """
    This class acts as a parser for .ins and .res SHELX files
    """
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

        # crystal site information
        self.crystal_sites = []
        self.mixed_site_numbers = []
        self.q_peaks = []

        self.r1 = 0.0
        self.suggested_weight_vals = ""
        self.read(filetxt)

    def read(self, filetxt):
        """
        Do actual parsing
        :param filetxt: name of file to read in
        :return the SHELXTL file
        """
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
        header_keys = ["TITL","CELL","ZERR","LATT","SYMM","SYMM","SYMM","SYMM","SYMM","SFAC","UNIT","TEMP","SIZE"]
        while True:
            if re.match("^\s*$", line) is None and all([key not in line for key in header_keys]):
                break
            line = lines[line_idx]
            self.extra_text[self.extra_text_section] += line + "\n"
            line_idx += 1

        line_idx += 1

        # Command/parameter section break on starting element key
        starting_element_keys = ["{}".format(el) for el in self.elements]
        while True:
            line = lines[line_idx]
            if re.match("^\s*$", line) is None:
                split = line.split()
                key = split[0]
                if key == "FVAR":
                    self.fvar_vals = [float(x) for x in split[1:]]
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


    def get_crystal_sites_text(self):
        return "\n".join([" ".join(cs.write_line()) for cs in self.crystal_sites])

    def get_ins_text(self):
        """
        Given a SHELXFile object, it returns the string that can be used to write it out
        :return: string of ins file
        """
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
        res += self.get_crystal_sites_text() + "\n"
        res += self.extra_text[2]
        return res

    # various editing methods ...
    def add_command(self, cmd, values=None):
        """
        Add a command to the ins file
        :param cmd: command to add
        """
        if cmd not in map(lambda tup: tup[0], self.commands):
            self.commands.append((cmd, values))

    def remove_command(self, cmd):
        """
        Remove command
        :param cmd: command to remove
        """
        self.commands = [(k, v) for k, v in self.commands if k != cmd]


    def add_anisotropy(self):
        self.add_command("ANIS")

    def remove_anisotropy(self):
        self.remove_command("ANIS")

    def add_exti(self):
        self.add_command("EXTI")

    def remove_exti(self):
        self.remove_command("EXTI")

    def get_site_stoichiometry(self, crystal_site):
        prefix = 1
        if crystal_site.occupancy_prefix < 0:
            prefix = 1.0 - self.fvar_vals[-1 * crystal_site.occupancy_prefix - 1]
        if crystal_site.occupancy_prefix > 1:
            prefix = self.fvar_vals[crystal_site.occupancy_prefix - 1]
        return crystal_site.occupancy * prefix

    def get_analytic_formula(self):
        formula = ""
        for site in self.crystal_sites:
            el = site.element
            stoich = self.get_site_stoichiometry(site)
            formula += el + stoich
        return formula

    def change_element(self, site_index, element_index):
        """
        Changes the element at a given site
        :param site_index: Index (starting at 0) of site
        :param element_index: Index (starting at 1) of element
        """
        site = self.get_crystal_sites_by_number(site_index + 1)[0]
        self.crystal_sites[site].element = element_index
        self.crystal_sites[site].el_string = self.elements[element_index - 1]
        self.crystal_sites[site].site_index = str(site_index+1)

    def move_q_to_crystal(self):
        """
        Move the top q peak to a crystal site
        """
        self.q_peaks[0].site_number = int(np.max(np.asarray([cs.site_number for cs in self.crystal_sites])) + 1.0)
        self.crystal_sites.append(self.q_peaks[0])
        del self.q_peaks[0]

    def move_crystal_to_q(self, site_idx=-1):
        """
        Move the crystal site (at site_idx) to a q peak
        :param site_idx: Index of crystal site to move
        :return:
        """
        site_idx = self.get_crystal_sites_by_number(site_idx)[0]
        self.q_peaks.insert(0, self.crystal_sites[site_idx])
        del self.crystal_sites[site_idx]

    def remove_sites_by_number(self, site_numbers):
        """
        Remove several crystal sites at once
        :param site_numbers: List of site numbers to remove
        :return:
        """
        self.crystal_sites = [site for site in self.crystal_sites if site.site_number not in site_numbers]

    def renumber_sites(self):
        """
        Renumbers the sites to consecutive integers
        :return:
        """
        renumbered = 0
        cur_site_number = None
        for site in self.crystal_sites:
            # Indexed from 1, so the first iteration will trigger this
            if cur_site_number != site.site_number:
                cur_site_number = site.site_number
                renumbered += 1
            site.site_number = renumbered

    def add_variable_occupancy(self, site_index):
        """
        If the site is currently fixed (with prefix=1), then unfix it by changing that prefix
        and adding an extra term to fvar_vals
        :param site_index: Index of site where it will unfix the occupancy
        """
        site = self.crystal_sites[site_index]
        if site.occupancy_prefix == 1:
            self.fvar_vals.append(0.5)
            site.occupancy_prefix = len(self.fvar_vals)

    def add_site_mixing(self, site_number, mixing_element_indices):
        """
        Adds multiple occupancy of a given crystal site.  Right now, only handles 2-element site mixing.
        :param site_number: Crystal site index
        :param mixing_element_indices: Indices of elements mixed at that site
        :return:
        """
        assert len(mixing_element_indices) <= 2, "Error: Can only handle mixing between 2 elements"
        site_indices = self.get_crystal_sites_by_number(site_number)
        self.commands.append(("EXYZ", [self.elements[i] + str(site_number) for i in mixing_element_indices]))
        self.commands.append(("EADP", [self.elements[i] + str(site_number) for i in mixing_element_indices]))

        # If site occupancy is not already variable, add a term to fvar:
        if self.crystal_sites[site_indices[0]].occupancy_prefix == 1:
            self.fvar_vals.append(0.5)

        mixed_sites = []
        for i, element_idx in enumerate(mixing_element_indices):
            new_site = copy.deepcopy(self.crystal_sites[site_indices[0]])
            new_site.el_string = self.elements[element_idx]
            new_site.site_number = str(site_number)
            new_site.element = element_idx + 1  # Because elements are 1-indexed
            if i == 0:
                new_site.occupancy_prefix = len(self.fvar_vals)
            else:
                new_site.occupancy_prefix = -len(self.fvar_vals)
            mixed_sites.append(new_site)

        self.remove_sites_by_number([site_number])  # Remove original crystal site
        self.crystal_sites.extend(mixed_sites)  # Replace with new mixed sites
        self.crystal_sites.sort(key=lambda site: site.site_number)
        self.mixed_site_numbers.append(site_number)

    def get_crystal_sites_by_number(self, number):
        """
        Returns the index in the crystal_site list corresponding to site_number=number
        In case of mixed sites, this list could have length > 1
        :param number: number of site we're looking for
        :return:
        """
        return [i for i, site in enumerate(self.crystal_sites) if site.site_number == number]

def main():
    test_file = "/Users/julialing/Documents/DataScience/crystal_refinement/4-2-1-4_single_crystal/Example_from_slides/7.ins"

    with open(test_file) as f:
        text = f.read()
        print text
        print "\n\n" + "~"*50 + "\n\n"
        file_obj = SHELXFile(text)

        print file_obj.write_ins()


if __name__ == "__main__":
    main()