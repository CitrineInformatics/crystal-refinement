import itertools
from pymatgen import Element
from pymatgen.structure_prediction.substitution_probability import SubstitutionProbability
from citrination_client import CitrinationClient


class OptimizerCache:
    """
    Class to hold persistent information for the optimizer.
    1) Bond length cache
    2) Mixing pairs
    """

    def __init__(self, shelx_file, bond_lengths=None, mixing_pairs=None, use_ml_model=False, api_key=None):
        # Chemistry information
        self.element_list = shelx_file.elements
        self.nominal_formula = shelx_file.get_nominal_formula()

        # Bond length information
        if use_ml_model:
            self.ml_model = CitrinationClient(api_key)
        else:
            self.ml_model = None

        self.bond_lengths = {}

        # Populate bond length cache
        for el1 in self.element_list:
            for el2 in self.element_list:
                self.get_bond_length(el1, el2)

        # Replace with any user-provided bond lengths
        if bond_lengths is not None:
            for el1, el2, bond_length in bond_lengths:
                self.bond_lengths[self._get_bond_key(el1, el2)] = bond_length

        # Mixing pairs information
        if mixing_pairs is None:
            self.mixing_pairs = self.get_mixing_pairs(probability_threshold=2E-4)
        else:
            self.mixing_pairs = []
            for el1, el2 in mixing_pairs:
                try:
                    self.mixing_pairs.append([shelx_file.get_element_by_name(el1),
                                              shelx_file.get_element_by_name(el2)])
                except KeyError:
                    "User defined mixing elements {} and {} must be elements specified in the ins file".format(el1, el2)

    @staticmethod
    def _get_bond_key(el1, el2):
        """
        Get the cache key for a pair of elements
        :param el1: element 1
        :param el2: element 2
        :return: cache key
        """
        return "-".join(sorted([el1, el2]))

    def get_bond_length(self, el1, el2):
        """
        Get the bond length as a function of the two elements and the formula. This
        :param el1: Element 1
        :param el2: Element 2
        :param formula: composition
        :return: predicted bond length
        """

        bond_key = self._get_bond_key(el1, el2)
        # Check the cache first
        if bond_key in self.bond_lengths:
            return self.bond_lengths[bond_key]

        bond_length = None

        if self.ml_model is not None:
            bond_length = self.get_predicted_bond_length(el1, el2)
        if bond_length is None:
            bond_length = self.get_naive_bond_length(el1, el2)

        self.bond_lengths[bond_key] = bond_length

        return bond_length

    @staticmethod
    def get_naive_bond_length(el1, el2):
        """
        Use pymatgen to get approximate bond length as the sum of atomic radii
        :param el1: Element 1
        :param el2: Element 2
        :return: estimated bond length
        """
        return Element(el1).atomic_radius + Element(el2).atomic_radius

    def get_predicted_bond_length(self, el1, el2):
        """
        Get the bond length prediction from the ML model as a function of the two elements and the formula. Since this
        prediction should be invariant to element order, it takes the average prediction from the two possible orders.
        If the uncertainty is too high, or there is an exception, this returns None
        :param el1: Element 1
        :param el2: Element 2
        :return: predicted bond length
        """

        formula = self.nominal_formula
        view_id = "680"

        try:
            candidate = [{"Element 1": el1, "Element 2": el2, "formula": formula},
                         {"Element 1": el2, "Element 2": el1, "formula": formula}]

            results = [x["Bond length"] for x in self.ml_model.predict(view_id, candidate)["candidates"]]

            # Only use the prediction if the uncertainty is low enough
            if results[0][1] < 0.5 or results[1][1] < 0.5:
                return (results[0][0] + results[1][0]) / 2

        except Exception:
            pass
        return None

    def get_shortest_bond(self):
        """
        Get the length of the shortest possible bond given the list of possible elements
        :return: length of the shortest possible bond
        """

        return sorted(self.bond_lengths.values())[0]

    def get_mixing_pairs(self, probability_threshold):
        """
        Get a sorted list of valid mixing pairs based on observed substitution probabilities.
        :param probability_threshold: probability threshold for designating a mixing pair as valid
        :return: sorted list of valid mixing pairs
        """
        pairs = []

        # For all elements in compound, calculate substitution probabilities
        # If substitution probability is > probability_threshold, then save it to pairs list.
        for e1, e2 in itertools.combinations(self.element_list, 2):
            sp = self.get_substitution_probability(e1.get_pymatgen_element(), e2.get_pymatgen_element())
            if sp > probability_threshold:
                pairs.append(([e1, e2], sp))
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
                total += sp.prob(self.get_ionic_specie(el1, o1), self.get_ionic_specie(el2, o2))
        return total

    @staticmethod
    def get_ionic_specie(el, ox):
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

    def get_report(self):
        """
        Return a report of the mixing pairs and bond lengths used, and whether they were calculated based on atomic
        radii or predicted by the ML model.
        This currently doesn't specify whether mixing pairs or bond lengths were provided by the user.
        :return: report as a string
        """
        report = ""
        report += "Mixing pairs considered: {}\n\n".format(
            ", ".join(["({}, {})".format(e1.get_name(), e2.get_name()) for e1, e2 in self.mixing_pairs]))

        if self.ml_model is None:
            report += "Bond lengths as calculated based on atomic radii:\n"
        else:
            report += "Bond lengths as predicted by Citrination machine learning model:\n"

        for bond_key, bond_length in sorted(self.bond_lengths.items(), key=lambda tup: tup[1]):
            report += "Bond: {}, length: {:.3f} ang\n".format(bond_key, bond_length)
        report += "\n"
        return report
