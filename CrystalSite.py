import numpy as np


class CrystalSite:
    def __init__(self, line_list):
        self.name = None
        self.position = None
        self.occupancy_prefix = None
        self.occupancy = None
        self.displacement = None
        self.electron_density = None
        self.anisotropy = None
        self.element = None
        self.read_line(line_list)

    def read_line(self, line_list):
        """
        Takes in a list of strings from the input file and sets the data members

        :param line_list: List of strings from the input file describing a crystal site
        """
        #Handle Q peaks:
        if (line_list[0][0] == "Q"):
            self.name = line_list[0]
            self.element = int(line_list[1])
            self.position = np.asarray([float(line_list[2]), float(line_list[3]), float(line_list[4])])
            self.occupancy_prefix = int(line_list[5][0])
            self.occupancy = float(line_list[5][1:])
            self.displacement = float(line_list[6])
            self.electron_density = float(line_list[7])

        # Handles assigned sites:
        else:
            self.name = line_list[0]
            self.element = int(line_list[1])
            self.position = np.asarray([float(line_list[2]), float(line_list[3]), float(line_list[4])])
            self.occupancy_prefix = int(line_list[5][0])
            self.occupancy = float(line_list[5][1:])
            self.displacement = float(line_list[6])
            if len(line_list) > 7:
                self.anisotropy = np.asarray([float(x) for x in line_list[7:]])

    def write_line(self):
        """
        Takes the crystal site data members and returns a list of strings which can then be written back out
        :return: line_list
        >>> cs = CrystalSite("SM2   4    0.153587    0.000000    0.432975    10.50000    0.00602".split())
        >>> cs.write_line()
        ['SM2', '4', '0.153587', '0.000000', '0.432975', '10.500000', '0.006020']

        """
        line_list = [self.name, str(self.element)]
        line_list += ['{:.6f}'.format(x) for x in self.position.tolist()]
        line_list += [str(self.occupancy_prefix) + '{:.6f}'.format(self.occupancy)]
        line_list += ['{:.6f}'.format(self.displacement)]
        if (self.name[0] == "Q"):
            line_list += ['{:.6f}'.format(self.electron_density)]
        if (self.anisotropy is not None):
            line_list += ['{:.6f}'.format(x) for x in self.anisotropy.tolist()]
        return line_list



    def set_element(self, element):
        self.element = element

    def set_position(self, position):
        self.position = position

    def get_position(self):
        return self.position

    def calc_min_distance_to_others(self, others):
        min_dist = np.infty
        for cs in others:
            dist = CrystalSite.calc_distance(self, cs)
            min_dist = np.minimum(min_dist, dist)
        return min_dist

    @staticmethod
    def calc_distance(site1, site2):
        """
        Calculate Euclidean distances between two sites
        :param site1: first site
        :param site2: second site
        :return:
        >>> cs1 = CrystalSite("SM2   4    0.153587    0.000000    0.432975    10.50000    0.00602".split())
        >>> cs2 = CrystalSite("IN3   3    0.000000   -0.500000    0.500000    20.25000    0.00749".split())
        >>> '{:.6f}'.format(CrystalSite.calc_distance(cs1, cs2))
        '0.527334'
        """
        assert(site1.position.shape == site2.position.shape)
        dist = np.sqrt(np.sum(np.square(site1.position - site2.position)))
        return dist

    @staticmethod
    def test_read_write():
        line_list = "SM1   4    0.086935    0.000000    0.881310    10.50000    0.00558".split()
        cs = CrystalSite(line_list)
        print cs.write_line()


if __name__ == "__main__":
    CrystalSite.test_read_write()

