import copy
import re

from itertools import groupby
from crystal_refinement.SHELXFile.SHELXElement import SHELXElement
from pymatgen import Composition

from crystal_refinement.SHELXFile.CrystalSite import CrystalSite


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
        # self.elements = []
        self._el_by_name = {}
        # self.formula_units = []

        # Stores commands as dict.  If command has no value, value is set to None
        self.commands = []
        self.fvar_vals = []

        # crystal site information
        self._crystal_sites = []
        self._mixed_site_indices = []
        self._crystal_sites_by_index = {}

        # self.crystal_sites = []
        # self.mixed_site_numbers = []
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
        # print filetxt

        # First text section, break at SFAC line
        while True:
            line = lines[line_idx]
            if "SFAC" in line:
                break
            self.extra_text[self.extra_text_section] += line + "\n"
            line_idx += 1

        # Parse composition info
        element_line = lines[line_idx].split()[1:]
        stoichiometry_line = lines[line_idx + 1].split()[1:]

        for i, (element, stoichiometry) in enumerate(zip(element_line, stoichiometry_line)):
            shelx_el = SHELXElement(element, stoichiometry, i + 1) # SHELX convention 1-indexes elements
            self.elements.append(shelx_el)
            self._el_by_name[element] = shelx_el
        # self.elements = lines[line_idx].split()[1:]
        # self.formula_units = lines[line_idx + 1].split()[1:]

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
        starting_element_keys = ["{}".format(el.get_name()) for el in self.elements]
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
                site_info = line.split()[:-1] + lines[line_idx + 1].split()
                line_idx += 1
            else:
                site_info = line.split()
            if site_info[0][0] == "Q":
                self.q_peaks.append(CrystalSite(site_info))
            else:
                self.add_site(CrystalSite(site_info))
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

    def get_element_string(self):
        res = "SFAC " + " ".join(map(lambda el: el.get_name(), self.elements)) + "\n"
        res += "UNIT " + " ".join(map(lambda el: el.nominal_stoichiometry, self.elements)) + "\n"
        return res

    def to_string(self):
        """
        Given a SHELXFile object, it returns the string that can be used to write it out
        :return: string of ins file
        """
        res = self.extra_text[0]
        res += self.get_element_string()
        res += self.extra_text[1] + "\n \n"
        for key, values in self.commands:
            res += key + "  "
            if values is not None:
                res += " ".join(values)
            res += "\n"
        res += "FVAR " + " ".join([str(x) for x in self.fvar_vals]) + "\n"
        res += self.crystal_sites_string() + "\n"
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

    def get_nominal_formula(self):
        formula = ""
        for el in self.elements:
            formula += "{}{}".format(el.get_name(True), el.nominal_stoichiometry)
        return Composition(formula).alphabetical_formula

    def get_analytic_formula(self):
        formula = ""
        for site in self.get_all_sites():
            el = site.el_string.capitalize()
            stoich = self.get_site_stoichiometry(site)
            formula += "{}{}".format(el, stoich)
        return Composition(formula).alphabetical_formula

    # def change_element(self, site_index, element_index):
    #     """
    #     Changes the element at a given site
    #     :param site_index: Index (starting at 0) of site
    #     :param element_index: Index (starting at 1) of element
    #     """
    #     site = self.crystal_sites_list.get_site_by_index(site_index)[0]
    #     site.element = element_index
    #     site.el_string = self.elements[element_index - 1]

    def move_q_to_crystal(self):
        """
        Move the top q peak to a crystal site
        """
        crystal_site = self.q_peaks[0]
        crystal_site.site_number = int(self.get_n_sites())
        self.add_site(crystal_site)
        self.q_peaks.remove(crystal_site)

    def move_crystal_to_q(self, site_idx=-1):
        """
        Move the crystal site (at site_idx) to a q peak
        :param site_idx: Index of crystal site to move
        :return:
        """

        if site_idx == -1:
            site_idx = self.get_n_sites() - 1
        sites = self.get_sites_by_index(site_idx)
        for site in sites:
            self.q_peaks.insert(0, site)
            self.remove_site(site)

    # def remove_sites_by_number(self, site_numbers):
    #     """
    #     Remove several crystal sites at once
    #     :param site_numbers: List of site numbers to remove
    #     :return:
    #     """
    #     self.crystal_sites = [site for site in self.crystal_sites if site.site_number not in site_numbers]

    # def renumber_sites(self):
    #     """
    #     Renumbers the sites to consecutive integers
    #     :return:
    #     """
    #     renumbered = 0
    #     cur_site_number = None
    #     for site in self.crystal_sites:
    #         # Indexed from 1, so the first iteration will trigger this
    #         if cur_site_number != site.site_number:
    #             cur_site_number = site.site_number
    #             renumbered += 1
    #         site.site_number = renumbered

    def add_variable_occupancy(self, site):
        """
        If the site is currently fixed (with prefix=1), then unfix it by changing that prefix
        and adding an extra term to fvar_vals
        :param site_index: Index of site where it will unfix the occupancy
        """

        if site.occupancy_prefix == 1:
            self.fvar_vals.append(0.5)
            site.occupancy_prefix = len(self.fvar_vals)

    def add_equal_site_mixing(self, site_index, mixing_elements):
        """
        Adds multiple occupancy of a given crystal site. Right now, only handles 2-element site mixing.
        Enforces equal 50:50 mixing between the elements.
        :param site_number: Crystal site index
        :param mixing_element_indices: Indices of elements mixed at that site
        :return:
        """
        assert len(mixing_elements) == 2, "Error: Can only handle mixing between 2 elements"
        sites = self.get_sites_by_index(site_index)
        assert len(sites) == 1, "Error: Site is already mixed"
        site = sites[0]
        # self.commands.append(("EXYZ", [self.elements[i] + str(site_index) for i in mixing_element_indices]))
        # self.commands.append(("EADP", [self.elements[i] + str(site_index) for i in mixing_element_indices]))
        self.commands.append(("EXYZ", [el.get_name() + str(site_index) for el in mixing_elements]))
        self.commands.append(("EADP", [el.get_name() + str(site_index) for el in mixing_elements]))
        for element in mixing_elements:
            new_site = copy.deepcopy(site)
            new_site.switch_element(element)
            new_site.occupancy = site.occupancy / 2
            self.add_site(new_site)
        self.remove_site(site)  # Remove original crystal site

    def add_site_mixing_variable_occupancy(self, site_index):
        """
        Sets the mixing on a site to be a free variable. Right now, only handles 2-element site mixing.
        :param site_number: Crystal site index
        :return:
        """
        sites = self.get_sites_by_index(site_index)
        # print(self.to_string())
        # for site in self.get_all_sites():
        #     print(site.to_string())
        # print(self._crystal_sites_by_index)
        assert len(sites) > 1, "Site {} must be mixed".format(site_index)
        # print(site_indices)
        # for site in self.crystal_sites:
        #     print(site.get_name())
        # If site occupancy is not already variable, add a term to fvar:
        if sites[0].occupancy_prefix == 1:
            self.fvar_vals.append(0.5)
            sites[0].occupancy_prefix = len(self.fvar_vals)
            sites[1].occupancy_prefix = -len(self.fvar_vals)
        else:
            sites[1].occupancy_prefix = -len(sites[0].occupancy_prefix)

    # def add_site_mixing(self, site_number, mixing_element_indices):
    #     """
    #     Adds multiple occupancy of a given crystal site.  Right now, only handles 2-element site mixing.
    #     :param site_number: Crystal site index
    #     :param mixing_element_indices: Indices of elements mixed at that site
    #     :return:
    #     """
    #     assert len(mixing_element_indices) <= 2, "Error: Can only handle mixing between 2 elements"
    #     site_indices = self.get_crystal_sites_by_number(site_number)
    #     self.commands.append(("EXYZ", [self.elements[i] + str(site_number) for i in mixing_element_indices]))
    #     self.commands.append(("EADP", [self.elements[i] + str(site_number) for i in mixing_element_indices]))
    #
    #     # If site occupancy is not already variable, add a term to fvar:
    #     if self.crystal_sites[site_indices[0]].occupancy_prefix == 1:
    #         self.fvar_vals.append(0.5)
    #
    #     mixed_sites = []
    #     for i, element_idx in enumerate(mixing_element_indices):
    #         new_site = copy.deepcopy(self.crystal_sites[site_indices[0]])
    #         new_site.el_string = self.elements[element_idx]
    #         new_site.site_number = str(site_number)
    #         new_site.element = element_idx + 1  # Because elements are 1-indexed
    #         if i == 0:
    #             new_site.occupancy_prefix = len(self.fvar_vals)
    #         else:
    #             new_site.occupancy_prefix = -len(self.fvar_vals)
    #         mixed_sites.append(new_site)
    #
    #     self.remove_sites_by_number([site_number])  # Remove original crystal site
    #     self.crystal_sites.extend(mixed_sites)  # Replace with new mixed sites
    #     self.crystal_sites.sort(key=lambda site: site.site_number)
    #     self.mixed_site_numbers.append(site_number)

    # def get_crystal_sites_by_number(self, index):
    #     """
    #     Returns the index in the crystal_site list corresponding to site_number=number
    #     In case of mixed sites, this list could have length > 1
    #     :param number: number of site we're looking for
    #     :return:
    #     """
    #     return self.crystal_sites_list.get_site_by_index(index)

    def _reindex_sites(self):
        """
        Repopulates the crystal site and mixed site indices. Also renumbers the sites to consecutive integers.
        :return:
        """
        self._crystal_sites_by_index = {}
        sorted_sites = sorted(self._crystal_sites, key=lambda site:site.site_number)
        for index, (_, sites) in enumerate(groupby(sorted_sites, key=lambda site: site.site_number)):
            site_list = []
            for site in sites:
                site.site_number = index + 1
                site_list.append(site)
            self._crystal_sites_by_index[index + 1] = site_list

        self._mixed_site_indices = []
        for index in self._crystal_sites_by_index.keys():
            if len(self._crystal_sites_by_index[index]) > 1:
                self._mixed_site_indices.append(index)

    # def reindex_sites_test(self):
    #     """
    #     Repopulates the crystal site and mixed site indices. Also renumbers the sites to consecutive integers.
    #     :return:
    #     """
    #     self._crystal_sites_by_index = {}
    #     groups = []
    #     for index, group in groupby(self._crystal_sites, key=lambda site: site.site_number):
    #         groups.append((index, list(group)))
    #
    #     for index, (_, sites) in enumerate(sorted(groups, key=lambda group: group[0])):
    #         site_list = []
    #         for site in list(sites):
    #             site.site_number = index + 1
    #             site_list.append(site)
    #         self._crystal_sites_by_index[index + 1] = site_list
    #
    #     self._mixed_site_indices = []
    #     for index in self._crystal_sites_by_index.keys():
    #         if len(self._crystal_sites_by_index[index]) > 1:
    #             self._mixed_site_indices.append(index)

    def add_site(self, new_site):
        self._crystal_sites.append(new_site)
        self._reindex_sites()

    def remove_site(self, site_to_remove):
        self._crystal_sites.remove(site_to_remove)
        self._reindex_sites()

    def crystal_sites_string(self):
        sorted_sites = sorted(self._crystal_sites, key=lambda site: site.site_number)
        return "\n".join(map(lambda site: site.to_string(), sorted_sites))

    def get_all_sites(self):
        return self._crystal_sites

    def get_sites_by_index(self, index):
        return self._crystal_sites_by_index[index]

    def get_mixed_site_indices(self):
        return self._mixed_site_indices

    def get_n_sites(self):
        return len(self._crystal_sites_by_index)

    def is_mixed_site(self, site):
        return site.site_number in self._mixed_site_indices

    def get_unmixed_sites(self):
        return filter(lambda site: not self.is_mixed_site(site), self._crystal_sites)

    def get_mixed_sites(self):
        return filter(lambda site: self.is_mixed_site(site), self._crystal_sites)

    def get_element_by_name(self, name):
        return self._el_by_name[name.upper()]

    def get_element_by_index(self, index):
        return self.elements[index]

    def copy(self):
        return copy.deepcopy(self)