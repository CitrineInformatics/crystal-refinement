from __future__ import absolute_import
from pymatgen import Composition
import crystal_refinement.utils.bond_utils as bond_utils
import math

"""
Various chemistry-based scores to evaluate a SHELX .ins/.res file. Smaller scores are better.
TODO: Define bounds for the different scores?
"""


def get_missing_element_penalty(shelx_file, cache):
    """
    Get the number of missing elements based on the nominal formula
    Bounds for this score are 0 to n_elements - 1.
    :param shelx_file: SHELX file
    :param cache: OptimizerCache object
    :return: missing elements score
    """
    nominal_elements = set(cache.element_names_list)
    result_elements = map(lambda site: site.get_element(True), shelx_file.get_all_sites())
    return len(nominal_elements.difference(result_elements))


def get_stoichiometry_score(shelx_file, cache):
    """
    Calculate a score for the given shelx file based on how well the composition agrees with the nominal formula
    Bounds for this score are 0 to 2.0 (maximum difference between two vectors of length 1).
    :param shelx_file: SHELX file
    :param cache: OptimizerCache object
    :return: stoichiometry agreement score
    """

    analytic_formula = Composition(shelx_file.get_analytic_formula())
    score = 0.0

    for el in cache.element_list:
        el_name = el.get_name(True)
        diff = abs(cache.nominal_formula.get_atomic_fraction(el_name) - analytic_formula.get_atomic_fraction(el_name))
        if diff > 0.05:
            score += diff
    return score


def get_compound_bond_score(bonds, shelx_file, cache):
    """
    Score a compound based on the bond lengths in the given SHELX file.
    Each site is given a score based on its nearest neighbor bonds, and the overall score is the stoichiometric-weighted
    average of these site scores.
    Bounds for this score are technically 0 to the longest nearest neighbor bond length.
    :param shelx_file: SHELX file
    :param cache: OptimizerCache object
    :return: bond score
    """

    site_bond_scores = bond_utils.get_site_bond_scores(bonds, cache, shelx_file, n_bonds=1)
    total_stoich = 0
    stoich_weighted_score = 0

    for site_name, score in site_bond_scores:
        stoich = 0
        for site in shelx_file.get_all_sites():
            if site.get_name().capitalize() == site_name:
                stoich = shelx_file.get_site_stoichiometry(site)
        total_stoich += stoich
        stoich_weighted_score += stoich * score

    return stoich_weighted_score / total_stoich


def get_negative_anisotropic_penalty(shelx_file):
    """
    Generate a penalty of 1.0 if any of the elements have a negative anisotropic displacement parameter in U11,
     U22 or U33
    Bounds for this score are 0 to 1.
    :param shelx_file: SHELx file
    :return: 1.0 for a file with negative ADPs, 0.0 otherwise
    """

    for site in shelx_file.get_all_sites():
        if site.anisotropy is not None:
            # Just check for negatives in U11, U22 and U33
            for adp in site.anisotropy[:3]:
                if adp < 0.0:
                    # If the ADP is negative, then return 1.0
                    return 1.0

    # If we haven't hit the early return, then all ADPs we checked were positive
    return 0.0


def get_overall_score(r1_score, bond_score, missing_elements_score, stoich_score, anisotropy_penalty, score_weighting):
    """
    Get the overall score for an optimizer result
    :param r1_score:
    :param bond_score:
    :param missing_elements_score:
    :param stoich_score:
    :param anisotropy_penalty:
    :param score_weighting:
    :return:
    """
    bond_score_basis = 0.1
    missing_elements_basis = 0.05
    stoich_score_basis = 0.1
    anisotropy_penalty_basis = 0.1

    # return math.pow(self.r1, self.score_weighting) + math.pow(self.bond_score, 1 - self.score_weighting)
    # return (1.0 + math.pow(self.r1 / r1_basis, self.score_weighting)) \
    #        * (1.0 + math.pow(self.bond_score / bond_score_basis, 1 - self.score_weighting))
    # try:
    return math.pow(r1_score, score_weighting) * \
           math.pow(bond_score + bond_score_basis, 1 - score_weighting) + \
           missing_elements_basis * missing_elements_score ** 2 + \
           stoich_score_basis * stoich_score + \
           anisotropy_penalty_basis * anisotropy_penalty
    # except TypeError:
    #     print(self.r1, self.score_weighting, self.bond_score, self.n_missing_elements)
    #     quit()
    # return self.r1 * self.bond_score
    # return (self.r1, len(self.res_file.mixed_site_numbers))