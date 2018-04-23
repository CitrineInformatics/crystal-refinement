from pymatgen import Composition
import crystal_refinement.utils.bond_utils as bond_utils

"""
Various chemistry-based scores to evaluate a SHELX .ins/.res file. Smaller scores are better.
"""


def get_missing_element_score(shelx_file, cache):
    """
    Get the number of missing elements based on the nominal formula
    :param shelx_file: SHELX file
    :param cache: OptimizerCache object
    :return: missing elements score
    """
    nominal_elements = cache.element_list
    result_elements = map(lambda site: site.el_string, shelx_file.get_all_sites())
    return len(nominal_elements.difference(result_elements))


def get_stoichiometry_score(shelx_file, cache):
    """
    Calculate a score for the given shelx file based on how well the composition agrees with the nominal formula
    :param shelx_file: SHELX file
    :param cache: OptimizerCache object
    :return: stoichiometry agreement score
    """

    analytic_formula = Composition(shelx_file.get_analytic_formula())
    score = 0.0
    # Is this the right way to iterate through?
    for el in cache.element_list:
        diff = abs(cache.nominal_formula.get_atomic_fraction(el) - analytic_formula.get_atomic_fraction(el))
        if diff > 0.05:
            score += diff
    return score


def get_compound_bond_score(bonds, shelx_file, cache):
    """
    Score a compound based on the bond lengths in the given SHELX file.
    Each site is given a score based on its nearest neighbor bonds, and the overall score is the stoichiometric-weighted
    average of these site scores.
    :param shelx_file: SHELX file
    :param cache: OptimizerCache object
    :return: bond score
    """

    site_bond_scores = bond_utils.get_site_bond_scores(bonds, cache, n_bonds=1)
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