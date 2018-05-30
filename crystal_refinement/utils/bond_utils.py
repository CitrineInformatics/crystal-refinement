from __future__ import absolute_import
import re, copy
from collections import defaultdict
from pymatgen.io.cif import CifParser


class Bond:
    def __init__(self, el1, el2, length):
        """
        :param el1: Element 1 as a SHELX site id
        :param el2: Element 2 as a SHELX site id
        :param length: length in angstroms
        """
        self.el1 = el1
        self.el2 = el2
        self.length = float(length)

    def get_normalized_element_names(self):
        """
        :return: tuple of the normalized element names
        """
        return [re.sub('\d', "", self.el1), re.sub('\d', "", self.el2)]

    def __cmp__(self, other):
        """
        Compare based on the length
        """
        if self.length < other.length:
            return -1
        if self.length == other.length:
            return 0
        if self.length > other.length:
            return 1

    def __str__(self):
        return "{}-{}: {}".format(self.el1, self.el2, self.length)


def get_site_bond_scores(bonds, cache, shelx_file, n_bonds=4):
    """
    Get a score for each site based on bond lengths.
    :param bonds: List of bond tuples (atom1, atom2, distance) for which to calculate priority
    :param cache: OptimizerCache object
    :param n_bonds: number of nearest neighbor bonds to incorporate in the score for each site
    :return: List of (site number, site score) tuples sorted with decreasing score ("worst" first).
    """
    nn_bonds_by_site = get_nn_bonds_by_site(bonds, shelx_file)

    bond_by_atom = defaultdict(lambda: [])
    for bond in nn_bonds_by_site.values():
        bond_score = get_bond_score(bond, nn_bonds_by_site, cache)
        bond_by_atom[bond.el1].append(bond_score)
        bond_by_atom[bond.el2].append(bond_score)
    # Average over scores from n_bonds shortest bonds
    res = map(lambda tup: (tup[0], sum(sorted(tup[1])[:n_bonds])), bond_by_atom.items())

    # Sort by which bonds are the shortest compared to what we'd expect
    return sorted(res, key=lambda tup: tup[1])


def get_bond_score(bond, nn_bonds_by_site, cache):
    """
    Get a score for a single bond.
    :param bond: The bond (el1, el2, length)
    :param nn_bonds_by_site: All nearest neighbor bonds in the compound
    :param cache: OptimizerCache object
    :return: bond score. Lower is better. For bonds that are ordered differently than ideal, the score will be
    appended with the difference between the bond lengths (some notion of how "wrong" it is).
    """
    score = 0.0
    ideal_bond_length = cache.get_ideal_bond_length(*bond.get_normalized_element_names())

    for other_bond in nn_bonds_by_site.values():
        other_ideal_bond_length = cache.get_ideal_bond_length(*other_bond.get_normalized_element_names())

        # Get the ideal and actual ordering for the bond lengths
        ideal_comp = ideal_bond_length < other_ideal_bond_length
        actual_comp = bond < other_bond

        # If the orderings are different, or the two bonds contain the same pair of elements
        if ideal_comp != actual_comp or ideal_bond_length == other_ideal_bond_length:
            # Difference between the two actual bond lengths
            diff = abs(bond.length - other_bond.length)
            ideal_diff = abs(ideal_bond_length - other_ideal_bond_length)
            if ideal_bond_length == other_ideal_bond_length:
                bond_score = (pow(max(diff - 0.1, 0.0), 2) + pow(max(ideal_diff - 0.1, 0.0), 2)) / 2
            else:
                bond_score = (pow(diff, 2) + pow(ideal_diff, 2)) / 2
            score += bond_score
    return score


def get_nn_bonds_by_site(bonds, shelx_file):
    """
    Get the first nearest neighbor bond associated with each site
    :param bonds: All bonds in the compound
    :param cache: OptimizerCache object
    :return: dictionary where the key is the site name and the value is the nearest neighbor bond for that site
    """

    sorted_bonds = sorted(bonds)

    nn_bonds_by_site = {}

    for site in shelx_file.get_all_sites():
        site_name = "{}{}".format(site.el_string, site.site_number)
        for bond in sorted_bonds:
            if site_name == bond.el1.upper() or site_name == bond.el2.upper():
                nn_bonds_by_site[site_name] = bond
                break
    return nn_bonds_by_site


def get_bonds(driver, shelx_file):
    """
    Get the bonds defined in the given shelx file
    :param driver: SHELXDriver object
    :param shelx_file: SHELXFile object
    :return: List of Bond objects
    """
    ins_file = copy.deepcopy(shelx_file)
    ins_file.add_command("ACTA")
    ins_file.remove_command("L.S.")
    ins_file.add_command("L.S.", ["1"])
    res = driver.run_SHELXTL(ins_file)
    if res is None:
        return []

    with open(driver.cif_file) as f:
        file_txt = f.read()
        cif_file = CifParser.from_string(file_txt)

    cif_dict = cif_file.as_dict().values()[0]
    bond_tuples = zip(cif_dict["_geom_bond_atom_site_label_1"], cif_dict["_geom_bond_atom_site_label_2"],
                      [float(x.replace("(", "").replace(")", "")) for x in cif_dict["_geom_bond_distance"]])

    return map(lambda tup: Bond(tup[0], tup[1], tup[2]), bond_tuples)
