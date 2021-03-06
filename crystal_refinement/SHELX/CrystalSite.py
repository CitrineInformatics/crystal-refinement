from __future__ import absolute_import
import numpy as np
import re


class CrystalSite:
    """
    Define class to hold information about a single crystal site or Q peak
    """
    def __init__(self, line_list):
        """
        Given a list of strings, assign values to the name, site_number, etc of the site

        :param line_list:
        """
        self.el_string = None # This would be FE if the name is FE1
        self.site_number = None   # This would be 1 if the name is FE1
        self.position = None  # This will be a 1X3 numpy array with the x, y, z coordinates
        self.occupancy_prefix = None  # This is a 1 in the case of fixed occupancy
        self.occupancy = None  # This number depends on the site symmetry.  Often it is 0.5
        self.displacement = None  # The displacement of the atom
        self.electron_density = None  # The calculated electron density
        self.anisotropy = None  # If anisotropy has been turned on, this records the anisotropy coefficients
        self.element = None  # This is the integer that refers to the element type in this location
        self.read_line(line_list)

    def read_line(self, line_list):
        """
        Takes in a list of strings from the input file and sets the data members

        :param line_list: List of strings from the input file describing a crystal site
        """
        #  Handle Q peaks:
        if line_list[0][0] == "Q":
            self.el_string = "Q"
            self.site_number = int(re.search('\d+', line_list[0]).group(0))
            self.element = int(line_list[1])
            self.position = np.asarray([float(line_list[2]), float(line_list[3]), float(line_list[4])])
            self.occupancy_prefix = int(line_list[5][0])
            self.occupancy = float(line_list[5][1:])
            self.displacement = float(line_list[6])
            self.electron_density = float(line_list[7])

        # Handles assigned sites:
        else:
            self.el_string = re.sub("\d+", "", line_list[0])
            self.site_number = int(re.search('\d+', line_list[0]).group(0))
            self.element = int(line_list[1])
            self.position = np.asarray([float(line_list[2]), float(line_list[3]), float(line_list[4])])
            self.occupancy_prefix = int(line_list[5][:(line_list[5].index(".") - 1)])
            self.occupancy = float(line_list[5][(line_list[5].index(".") - 1):])
            self.displacement = float(line_list[6])
            if len(line_list) > 7:
                self.anisotropy = np.asarray([float(x) for x in line_list[7:]])

    def to_string(self):
        """
        Takes the crystal site data members and returns a list of strings which can then be written back out
        :return: line_list
        >>> cs = CrystalSite("SM2   4    0.153587    0.000000    0.432975    10.50000    0.00602".split())
        >>> cs.to_string()
        ['SM2', '4', '0.153587', '0.000000', '0.432975', '10.500000', '0.006020']

        """
        line_list = [self.el_string + str(self.site_number), str(self.element)]
        line_list += ['{:.6f}'.format(x) for x in self.position.tolist()]
        line_list += [str(self.occupancy_prefix) + '{:.6f}'.format(self.occupancy)]
        line_list += ['{:.6f}'.format(self.displacement)]
        if (self.el_string == "Q"):
            line_list += ['{:.6f}'.format(self.electron_density)]
        if (self.anisotropy is not None):
            line_list += ['=\n'] + ['{:.6f}'.format(x) for x in self.anisotropy.tolist()]
        return " ".join(line_list)

    def get_name(self):
        """
        Get the string representation of the site name, which is the element string + site number
        :return: site name as a string
        """
        return self.el_string + str(self.site_number)

    def get_element(self, capitalized=False):
        """
        Get the element name for this crystal site
        :param capitalized: Whether to capitalize conventionally (Helium would be He)
        :return: element name
        """
        if capitalized:
            return self.el_string.capitalize()
        else:
            return self.el_string

    def set_element(self, element):
        """
        Set the element based on the index
        :param element: element index to switch to
        """
        self.element = element

    def get_position(self):
        """
        Get the position tuple
        :return: the position of the site as a tuple
        """
        return self.position

    def set_position(self, position):
        """
        Set the position tuple
        """
        self.position = position

    def switch_element(self, shelx_element):
        """
        Get the position tuple
        :return: the position of the site as a tuple
        """
        self.element = shelx_element.index
        self.el_string = shelx_element.name




