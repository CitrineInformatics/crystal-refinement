from pymatgen.structure_prediction.substitution_probability import SubstitutionProbability
from pymatgen.core.composition import Element
import re
from collections import defaultdict


def get_specie(el, ox):
    if ox > 0:
        return el.symbol + str(ox) + "+"
    else:
        return el.symbol + str(abs(ox)) + "-"


def get_substitution_probability(el1, el2):
    sp = SubstitutionProbability()
    total = 0
    for o1 in el1.common_oxidation_states:
        for o2 in el2.common_oxidation_states:
            total += sp.prob(get_specie(el1, o1), get_specie(el2, o2))
    return total


def site_mixing_priority(bonds, use_ml=False):
    site_scores = get_site_scores(bonds, use_ml)
    return map(lambda tup: int(re.search("\d+", tup[0]).group(0)), site_scores)


def get_site_scores(bonds, use_ml=False):
    bond_by_atom = defaultdict(lambda: [])
    for a1, a2, distance in bonds:
        ideal_bond_length = get_ideal_bond_length(a1, a2, use_ml)
        bond_by_atom[a1].append(distance - ideal_bond_length)
        bond_by_atom[a2].append(distance - ideal_bond_length)
    return sorted(map(lambda tup: (tup[0], sum(sorted(tup[1])[:4])), bond_by_atom.items()), key=lambda tup: tup[1])


def get_ideal_bond_length(specie_name1, specie_name2, use_ml=False):
    if use_ml:
        pass
    else:
        el1 = Element(re.sub('\d', "", specie_name1))
        el2 = Element(re.sub('\d', "", specie_name2))
        return el1.atomic_radius + el2.atomic_radius
