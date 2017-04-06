from pymatgen.structure_prediction.substitution_probability import SubstitutionProbability
from pymatgen.core.composition import Element
import re
from collections import defaultdict


def get_specie(el, ox):
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


def get_substitution_probability(el1, el2):
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
            total += sp.prob(get_specie(el1, o1), get_specie(el2, o2))
    return total


def get_bond_score(bond, use_ml_model=False):
    """
    Scores the likelihood of a bond based on its deviation from an ideal bond length.
    :param bond: The bond to score
    :param use_ml_model: Whether to use an ml model to calculate the ideal bond length.
    TODO: The ML model needs access to the .ins file to get the chemical system
    TODO: Might want to refactor this to use uncertainty estimates to determine the likelihood of the given bond length
    If not, the ideal bond length is calculated as the sum of atomic radii.
    TODO: This should be covalent radii instead
    :return: bond length - ideal bond length
    """
    ideal_bond_length = get_ideal_bond_length(bond[0], bond[1], use_ml_model)
    # Calculate amount which bonds are shorter than they are supposed to be
    return bond[2] - ideal_bond_length


def site_mixing_priority(bonds, use_ml_model=False):
    """
    Determines priority for adding in site mixing.  Finds most problematic sites based on whether bonds
    are shorter than expected.
    TODO: What if bonds are longer than expected?
    :param bonds: List of bond tuples (atom1, atom2, distance) for which to calculate priority
    :return: List of (site number, site score) tuples sorted with decreasing priority.
    """
    bond_by_atom = defaultdict(lambda: [])
    for bond in bonds:
        bond_score = get_bond_score(bond, use_ml_model)
        bond_by_atom[bond[0]].append(bond_score)
        bond_by_atom[bond[1]].append(bond_score)
    # Average over scores from 4 shortest bonds
    res = map(lambda tup: (tup[0], sum(sorted(tup[1])[:4])), bond_by_atom.items())
    # Sort by which bonds are the shortest compared to what we'd expect
    res = sorted(res, key=lambda tup: tup[1])
    # Get associated site index
    return map(lambda tup: (int(re.search("\d+", tup[0]).group(0)), tup[1]), res)


def get_ideal_bond_length(specie_name1, specie_name2, use_ml_model=False):
    """
    Scores the likelihood of a bond based on its deviation from an ideal bond length.
    :param specie_name1: One element in the bond
    :param specie_name2: The other element in the bond
    :param use_ml_model: Whether to use an ml model to calculate the ideal bond length.
    TODO: The ML model needs access to the .ins file to get the chemical system
    If not, the ideal bond length is calculated as the sum of atomic radii.
    TODO: This should be covalent radii instead
    :return: ideal bond length
    """
    if use_ml_model:
        pass
    else:
        # Use pymatgen to get approximate bond length = sum of atomic radii
        el1 = Element(re.sub('\d', "", specie_name1))
        el2 = Element(re.sub('\d', "", specie_name2))
        return el1.atomic_radius + el2.atomic_radius
