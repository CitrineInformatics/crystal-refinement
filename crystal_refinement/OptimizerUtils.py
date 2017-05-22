from pymatgen.structure_prediction.substitution_probability import SubstitutionProbability
from pymatgen.core.composition import Element
import re, copy, itertools, os
from collections import defaultdict
from pymatgen.io.cif import CifParser
from citrination_client import CitrinationClient

class OptimizerUtils:
    """
    Class for performing single crystal refinement
    """

    def __init__(self, original_ins, bond_lengths=None, mixing_pairs=None, use_ml_model=False):
        if use_ml_model:
            self.ml_model = CitrinationClient(os.environ["CITRINATION_API_KEY"])
        else:
            self.ml_model = None
        self.original_ins = original_ins
        if bond_lengths is None:
            self.bond_lengths = {}
        else:
            self.bond_lengths = {self.get_bond_key(el1, el2) for el1, el2, bond_length in bond_lengths}
        if mixing_pairs is None:
            self.mixing_pairs = self.get_mixing_pairs(probability_threshold=2E-4)
        else:
            self.mixing_pairs = mixing_pairs

    def get_mixing_pairs(self, probability_threshold):
        element_list = [Element(el.capitalize()) for el in self.original_ins.elements]
        pairs = []

        # For all elements in compound, calculate substitution probabilities
        # If substitution probability is > probability_threshold, then save it to pairs list.
        for i1, i2 in itertools.combinations(range(len(element_list)), 2):
            e1 = element_list[i1]
            e2 = element_list[i2]
            sp = self.get_substitution_probability(e1, e2)
            # print e1, e2, sp
            if sp > probability_threshold:
                pairs.append(([i1, i2], sp))

        # Sort pairs by substitution probability (largest to smallest)
        return [tup[0] for tup in sorted(pairs, key=lambda tup: -tup[1])]

    def get_substitution_probability(self, el1, el2):
        """
        Use pymatgen to calculate substitution probability for use in site mixing
        :param el1: name of first element
        :param el2: name of second element
        :return: total probability of substitution between two elements
        """
        sp = SubstitutionProbability()
        total = 0
        for o1 in el1.common_oxidation_states:
            for o2 in el2.common_oxidation_states:
                total += sp.prob(self.get_specie(el1, o1), self.get_specie(el2, o2))
        return total

    def get_shortest_bond(self, ins_file):
        return sorted([self.get_ideal_bond_length(el.capitalize(), el.capitalize()) for el in ins_file.elements])[0]

    def get_specie(self, el, ox):
        """
        Return string for ion
        :param el: Element abbreviation (eg. Fe, Au, ...)
        :param ox: Oxidation number (eg. +2)
        :return: eg. Fe2+
        """
        if ox > 0:
            return el.symbol + str(ox) + "+"
        else:
            return el.symbol + str(abs(ox)) + "-"


    def score_compound_bonds(self, bonds, shelx_file):
        """
        The smaller the better
        :param bonds:
        :param ml_model:
        :return:
        """
        site_bond_scores = self.get_site_bond_scores(bonds, 1)
        total_stoich = 0
        stoich_weighted_score = 0
        for site_name, score in site_bond_scores:
            stoich = 0
            for site in shelx_file.crystal_sites:
                if site.get_name().capitalize() == site_name:
                    stoich = shelx_file.get_site_stoichiometry(site)
            total_stoich += stoich
            stoich_weighted_score += stoich * score

        if len(site_bond_scores) < len(shelx_file.crystal_sites):
            stoich_weighted_score -= 50

        return -stoich_weighted_score / total_stoich

    def get_bond_score(self, bond):
        """
        Scores the likelihood of a bond based on its deviation from an ideal bond length.
        :param bond: The bond to score
        :param ml_model: Whether to use an ml model to calculate the ideal bond length.
        TODO: The ML model needs access to the .ins file to get the chemical system
        TODO: Might want to refactor this to use uncertainty estimates to determine the likelihood of the given bond length
        If not, the ideal bond length is calculated as the sum of atomic radii.
        TODO: This should be covalent radii instead
        :return: bond length - ideal bond length
        """
        ideal_bond_length = self.get_ideal_bond_length(bond[0], bond[1])
        # Calculate amount which bonds deviate from the expected bond length
        return -abs(bond[2] - ideal_bond_length)

    def get_bond_score2(self, bond):
        """
        Scores the likelihood of a bond based on its deviation from an ideal bond length.
        :param bond: The bond to score
        :param ml_model: Whether to use an ml model to calculate the ideal bond length.
        TODO: The ML model needs access to the .ins file to get the chemical system
        TODO: Might want to refactor this to use uncertainty estimates to determine the likelihood of the given bond length
        If not, the ideal bond length is calculated as the sum of atomic radii.
        TODO: This should be covalent radii instead
        :return: bond length - ideal bond length
        """
        ideal_bond_length = self.get_ideal_bond_length(bond[0], bond[1])
        # Calculate amount which bonds deviate from the expected bond length
        return -((bond[2] - ideal_bond_length) / 2) ** 2

    def get_site_bond_scores(self, bonds, n_bonds=4):
        """
            Determines priority for adding in site mixing.  Finds most problematic sites based on whether bonds
            are shorter than expected.
            TODO: What if bonds are longer than expected?
            :param bonds: List of bond tuples (atom1, atom2, distance) for which to calculate priority
            :return: List of (site number, site score) tuples sorted with decreasing priority.
            """
        bond_by_atom = defaultdict(lambda: [])
        for bond in bonds:
            bond_score = self.get_bond_score(bond)
            bond_by_atom[bond[0]].append(bond_score)
            bond_by_atom[bond[1]].append(bond_score)
        # Average over scores from 4 shortest bonds
        res = map(lambda tup: (tup[0], sum(sorted(tup[1])[:n_bonds])), bond_by_atom.items())
        # Sort by which bonds are the shortest compared to what we'd expect
        return sorted(res, key=lambda tup: tup[1])

    def site_mixing_priority(self, bonds, n_bonds=4):

        # Get associated site index
        return map(lambda tup: (int(re.search("\d+", tup[0]).group(0)), tup[1]), self.get_site_bond_scores(bonds, n_bonds))

    def get_bond_key(self, el1, el2):
        return ",".join(sorted([el1, el2]))

    def get_ideal_bond_length(self, specie_name1, specie_name2):
        """
        Scores the likelihood of a bond based on its deviation from an ideal bond length.
        :param specie_name1: One element in the bond
        :param specie_name2: The other element in the bond
        TODO: The ML model needs access to the .ins file to get the chemical system
        If not, the ideal bond length is calculated as the sum of atomic radii.
        TODO: This should be covalent radii instead
        :return: ideal bond length
        """

        bond_key = self.get_bond_key(re.sub('\d', "", specie_name1), re.sub('\d', "", specie_name2))
        if bond_key in self.bond_lengths:
            return self.bond_lengths[bond_key]

        if self.ml_model is not None:
            pass

        # Use pymatgen to get approximate bond length = sum of atomic radii
        el1 = Element(re.sub('\d', "", specie_name1))
        el2 = Element(re.sub('\d', "", specie_name2))
        return el1.atomic_radius + el2.atomic_radius


    def get_bonds(self, driver, ins_file):
        """
        Get the bonds defined in the given ins_file

        :param ins_file: SHELXFile object
        :return res: List of bond tuples (element 1, element 2, bond length)
        """
        ins_file = copy.deepcopy(ins_file)
        ins_file.add_command("ACTA")
        ins_file.remove_command("L.S.")
        ins_file.add_command("L.S.", ["1"])
        driver.run_SHELXTL(ins_file)

        with open(driver.cif_file) as f:
            cif_file = CifParser.from_string(f.read())

        cif_dict = cif_file.as_dict().values()[0]
        return zip(cif_dict["_geom_bond_atom_site_label_1"], cif_dict["_geom_bond_atom_site_label_2"],
                   [float(x.replace("(", "").replace(")", "")) for x in cif_dict["_geom_bond_distance"]])
