# from pymatgen.structure_prediction.substitution_probability import SubstitutionProbability
# from pymatgen import Element, Composition
# import re, copy, itertools, os
# from collections import defaultdict
# from pymatgen.io.cif import CifParser
# from citrination_client import CitrinationClient
# import warnings
#
# class OptimizerUtils:
#     """
#     Class for performing single crystal refinement
#     """
#
#     def __init__(self, shelx_file, bond_lengths=None, mixing_pairs=None, use_ml_model=False):
#         self.element_list = [str(Element(el.get_name(True))) for el in shelx_file.elements]
#         if use_ml_model:
#             try:
#                 self.ml_model = CitrinationClient(os.environ["CITRINATION_API_KEY"])
#             except KeyError:
#                 warnings.warn("To use the machine learning model, you need a Citrination account.")
#                 print("To use the machine learning bond length model, you need a free account from citrination.com. " +\
#                       "Please set the environment variable " + \
#                       "CITRINATION_API_KEY using your Citrination API key. ")
#                 print("Defaulting to simple atomic radii bond length model.")
#                 self.ml_model = None
#             self.prediction_cache = {}
#         else:
#             self.ml_model = None
#         self.bond_lengths = {}
#
#         for el1 in self.element_list:
#             for el2 in self.element_list:
#                 bond_key = self.get_bond_key(el1, el2)
#                 if bond_key not in self.bond_lengths.keys():
#                     self.bond_lengths[bond_key] = self.get_bond_length(el1, el2, shelx_file)
#
#         if bond_lengths is not None:
#             for el1, el2, bond_length in bond_lengths:
#                 self.bond_lengths[self.get_bond_key(el1, el2)] = bond_length
#
#         if mixing_pairs is None:
#             self.mixing_pairs = self.get_mixing_pairs(shelx_file, probability_threshold=2E-4)
#         else:
#             self.mixing_pairs = []
#             for el1, el2 in mixing_pairs:
#                 try:
#                     self.mixing_pairs.append([shelx_file.get_element_by_name(el1),
#                                               shelx_file.get_element_by_name(el2)])
#                 except KeyError:
#                     "User defined mixing elements {} and {} must be elements specified in the ins file".format(el1, el2)
#
#     def get_ml_prediction(self, el1, el2, shelx_file):
#         """
#         Get the bond length prediction from the ML model as a function of the two elements and the formula. Since this
#         prediction should be invariant to element order, it takes the average prediction from the two possible orders.
#         :param el1: Element 1
#         :param el2: Element 2
#         :param shelx_file: shelx file which specifies the formula.
#         :return:
#         """
#         formula = shelx_file.get_analytic_formula()
#         prediction_key = (self.get_bond_key(el1, el2), formula)
#         # Check the cache first
#         if prediction_key in self.prediction_cache:
#             return self.prediction_cache[prediction_key][0]
#         try:
#             candidate = [{"Element 1": el1, "Element 2": el2, "formula": formula},
#                          {"Element 1": el2, "Element 2": el1, "formula": formula}]
#
#             results = [x["Bond length"] for x in self.ml_model.predict("680", candidate)["candidates"]]
#             uncertainty = (results[0][1] + results[1][1]) / 2
#             # Only use the prediction if the uncertainty is low enough
#             if uncertainty < 0.5:
#                 bond_length_prediction = (results[0][0] + results[1][0]) / 2
#                 self.prediction_cache[prediction_key] = [bond_length_prediction, uncertainty]
#                 return self.prediction_cache[prediction_key][0]
#         except Exception:
#             pass
#         self.prediction_cache[prediction_key] = [float(Element(el1).atomic_radius + Element(el2).atomic_radius), 0.0]
#         return self.prediction_cache[prediction_key][0]
#
#     def get_report(self):
#         """
#         Return a report of the mixing pairs and bond lengths used, and whether they were calculated basd on atomic
#         radii or predicted by the ML model.
#         This currently doesn't specify whether mixing pairs or bond lengths were provided by the user.
#         :return: report as a string
#         """
#         report = ""
#         report += "Mixing pairs considered: {}\n\n".format(", ".join(["({}, {})".format(e1.get_name(), e2.get_name()) for e1, e2 in self.mixing_pairs]))
#
#         if self.ml_model is None:
#             report += "Bond lengths as calculated based on atomic radii:\n"
#         else:
#             report += "Bond lengths as predicted by Citrination machine learning model:\n"
#
#         for bond_key, bond_length in sorted(self.bond_lengths.items(), key=lambda tup: tup[1]):
#             report += "Bond: {}, length: {:.3f} ang\n".format(bond_key, bond_length)
#         report += "\n"
#         return report
#
#     def get_mixing_pairs(self, shelx_file, probability_threshold):
#         """
#         Get a sorted list of valid mixing pairs based on observed substitution probabilities.
#         :param shelx_file: SHELX file that specifies the list of potential elements.
#         :param probability_threshold: probability threshold for designating a mixing pair as valid
#         :return: sorted list of valid mixing pairs
#         """
#         pairs = []
#
#         # For all elements in compound, calculate substitution probabilities
#         # If substitution probability is > probability_threshold, then save it to pairs list.
#         for e1, e2 in itertools.combinations(shelx_file.elements, 2):
#             sp = self.get_substitution_probability(e1.get_pymatgen_element(), e2.get_pymatgen_element())
#             if sp > probability_threshold:
#                 pairs.append(([e1, e2], sp))
#         # Sort pairs by substitution probability (largest to smallest)
#         return [tup[0] for tup in sorted(pairs, key=lambda tup: -tup[1])]
#
#     def get_substitution_probability(self, el1, el2):
#         """
#         Use pymatgen to calculate substitution probability for use in site mixing
#         :param el1: name of first element
#         :param el2: name of second element
#         :return: total probability of substitution between two elements
#         """
#         sp = SubstitutionProbability()
#         total = 0
#         for o1 in el1.common_oxidation_states:
#             for o2 in el2.common_oxidation_states:
#                 total += sp.prob(self.get_specie(el1, o1), self.get_specie(el2, o2))
#         return total
#
#     def get_shortest_bond(self, shelx_file):
#         """
#         Get the length of the shortest possible bond given the list of elements
#         :param shelx_file: SHELX file that specifies the list of potential elements.
#         :return: length of the shortest possible bond
#         """
#         return sorted([self.get_ideal_bond_length(el.get_name(True), el.get_name(True), shelx_file) for el in shelx_file.elements])[0]
#
#     def get_specie(self, el, ox):
#         """
#         Return string for ion
#         :param el: Element abbreviation (eg. Fe, Au, ...)
#         :param ox: Oxidation number (eg. +2)
#         :return: eg. Fe2+
#         """
#         if ox > 0:
#             return el.symbol + str(ox) + "+"
#         else:
#             return el.symbol + str(abs(ox)) + "-"
#
#
#     def score_compound_bonds(self, bonds, shelx_file):
#         """
#         The smaller the better
#         :param bonds:
#         :param ml_model:
#         :return:
#         """
#         site_bond_scores = self.get_site_bond_scores(bonds, shelx_file, n_bonds=1)
#         total_stoich = 0
#         stoich_weighted_score = 0
#         for site_name, score in site_bond_scores:
#             stoich = 0
#             for site in shelx_file.crystal_sites:
#                 if site.get_name().capitalize() == site_name:
#                     stoich = shelx_file.get_site_stoichiometry(site)
#             total_stoich += stoich
#             stoich_weighted_score += stoich * score
#
#         if len(site_bond_scores) < len(shelx_file.crystal_sites):
#             stoich_weighted_score -= 50
#
#         return -stoich_weighted_score / total_stoich
#
#     def get_bond_score(self, bond, shelx_file):
#         """
#         Scores the likelihood of a bond based on its deviation from an ideal bond length.
#         :param bond: The bond to score
#         :param ml_model: Whether to use an ml model to calculate the ideal bond length.
#         TODO: The ML model needs access to the .ins file to get the chemical system
#         TODO: Might want to refactor this to use uncertainty estimates to determine the likelihood of the given bond length
#         If not, the ideal bond length is calculated as the sum of atomic radii.
#         TODO: This should be covalent radii instead
#         :return: bond length - ideal bond length
#         """
#         ideal_bond_length = self.get_ideal_bond_length(bond[0], bond[1], shelx_file)
#         # Calculate amount which bonds deviate from the expected bond length
#         return -abs(bond[2] - ideal_bond_length)
#
#     def get_bond_score2(self, bond, shelx_file):
#         """
#         Scores the likelihood of a bond based on its deviation from an ideal bond length.
#         :param bond: The bond to score
#         :param ml_model: Whether to use an ml model to calculate the ideal bond length.
#         :return: bond length - ideal bond length
#         """
#         ideal_bond_length = self.get_ideal_bond_length(bond[0], bond[1], shelx_file)
#         # Calculate amount which bonds deviate from the expected bond length
#         return -((bond[2] - ideal_bond_length) / 2) ** 2
#
#     def get_bond_score_relative_distance(self, bond, all_nn_bonds, verbose=False):
#         """
#         :param bond: The bond (el1, el2, length)
#         :param all_nn_bonds: All nearest neighbor bonds in the compound
#         :param shelx_file:
#         :return: bond score. Lower is better. For bonds that are ordered differently than ideal, the score will be
#         appended with the difference between the bond lengths (some notion of how "wrong" it is).
#         """
#         if verbose:
#             print("calculating bond score for {}".format(bond))
#         score = 0.0
#         ideal_bond_length = self.bond_lengths[self.get_bond_key(self.specie_to_el(bond[0]), self.specie_to_el(bond[1]))]
#         for other_bond in all_nn_bonds:
#             other_ideal_bond_length = self.bond_lengths[self.get_bond_key(self.specie_to_el(other_bond[0]), self.specie_to_el(other_bond[1]))]
#             ideal_comp = ideal_bond_length < other_ideal_bond_length
#             actual_comp = bond[2] < other_bond[2]
#             if ideal_comp != actual_comp or ideal_bond_length == other_ideal_bond_length:
#                 diff = abs(bond[2] - other_bond[2])
#                 ideal_diff = abs(ideal_bond_length - other_ideal_bond_length)
#                 if ideal_bond_length == other_ideal_bond_length:
#                     bond_score = (pow(max(diff - 0.1, 0.0), 2) + pow(max(ideal_diff - 0.1, 0.0), 2)) / 2
#                 else:
#                     bond_score = (pow(diff, 2) + pow(ideal_diff, 2)) / 2
#                 if verbose and bond_score > 0.0:
#                     print("Bond added to score: {}, {}, {}, {}".format(diff, bond_score, bond, other_bond))
#                 score += bond_score
#         if verbose:
#             print("final score for {} is {}".format(bond, score))
#         return score
#
#     def get_site_bond_scores_relative_distance(self, bonds, all_nn_bonds, n_bonds=4, verbose=False):
#         """
#         Determines priority for adding in site mixing.  Finds most problematic sites based on whether bonds
#         are shorter than expected.
#         TODO: What if bonds are longer than expected?
#         :param bonds: List of bond tuples (atom1, atom2, distance) for which to calculate priority
#         :return: List of (site number, site score) tuples sorted with decreasing priority.
#         """
#         bond_by_atom = defaultdict(lambda: [])
#         for bond in bonds:
#             bond_score = self.get_bond_score_relative_distance(bond, all_nn_bonds, verbose=verbose)
#             bond_by_atom[bond[0]].append(bond_score)
#             bond_by_atom[bond[1]].append(bond_score)
#         # Average over scores from n shortest bonds
#         res = map(lambda tup: (tup[0], sum(sorted(tup[1])[:n_bonds])), bond_by_atom.items())
#         if verbose:
#             print("scores after filtering to shortest bonds:")
#             for bond in res:
#                 print(bond)
#         # Sort by which bonds are the shortest compared to what we'd expect
#         return sorted(res, key=lambda tup: tup[1])
#
#     def get_nn_bonds_by_site(self, bonds, shelx_file):
#         sorted_bonds = sorted(bonds, key=lambda tup: tup[2])
#
#         nn_bonds_by_site = {}
#
#         # print(len(shelx_file.crystal_sites))
#         for site in shelx_file.get_all_sites():
#             site_name = "{}{}".format(site.el_string, site.site_number)
#             # print site_name
#             for bond in sorted_bonds:
#                 if site_name == bond[0].upper() or site_name == bond[0].upper():
#                     nn_bonds_by_site[site_name] = bond
#                     # print(bond)
#                     break
#         return nn_bonds_by_site
#
#     def score_compound_bonds_relative_distance(self, bonds, shelx_file, driver, verbose=False):
#         """
#         The smaller the better
#         :param bonds:
#         :param ml_model:
#         :return:
#         """
#
#         nn_bonds_by_site = self.get_nn_bonds_by_site(bonds, shelx_file)
#
#
#         # for bond in bonds:
#         #     print(bond)
#
#         # print("\n\n\n")
#
#         # print("NN bonds:")
#         # for site, bond in nn_bonds_by_site.items():
#         #     print(site, bond)
#         # quit()
#
#         if verbose:
#             print("bonds by site:")
#             for bond in nn_bonds_by_site.items():
#                 print(bond)
#
#         site_bond_scores = self.get_site_bond_scores_relative_distance(nn_bonds_by_site.values(), nn_bonds_by_site.values(), n_bonds=1, verbose=verbose)
#         total_stoich = 0
#         stoich_weighted_score = 0
#
#         for site_name, score in site_bond_scores:
#             if verbose:
#                 print(site_name, score)
#             stoich = 0
#             for site in shelx_file.get_all_sites():
#                 if site.get_name().capitalize() == site_name:
#                     stoich = shelx_file.get_site_stoichiometry(site)
#             total_stoich += stoich
#             stoich_weighted_score += stoich * score
#
#         # print("crystal sites:")
#         # for site in shelx_file.crystal_sites:
#         #     print(site.get_name().capitalize())
#         # print("bonds passed in:")
#         # for bond in bonds:
#         #     print(bond)
#         # if len(nn_bonds_by_site) == 0:
#             # print("bonds passed in:")
#             # for bond in bonds:
#             #     print(bond)
#             # print("bonds again:")
#             # bonds_again = self.get_bonds(driver, shelx_file)
#             # for bond in bonds_again:
#             #     print(bond)
#             # bond_by_atom = defaultdict(lambda: [])
#             # for bond in nn_bonds_by_site.values():
#             #     bond_score = self.get_bond_score_relative_distance(bond, nn_bonds_by_site.values())
#             #     bond_by_atom[bond[0]].append(bond_score)
#             #     bond_by_atom[bond[1]].append(bond_score)
#             # print(len(nn_bonds_by_site.values()))
#             # print(len(bond_by_atom))
#             # Average over scores from n shortest bonds
#             # res = map(lambda tup: (tup[0], sum(sorted(tup[1])[:1])), bond_by_atom.items())
#             # Sort by which bonds are the shortest compared to what we'd expect
#             # print len(shelx_file.crystal_sites)
#             # print(site_bond_scores)
#             # print(stoich_weighted_score, total_stoich)
#             # print("{} bonds".format(len(bonds)))
#             # print(shelx_file.get_ins_text())
#
#         # if len(site_bond_scores) < len(shelx_file.crystal_sites):
#             # stoich_weighted_score -= 50
#         # print len(shelx_file.crystal_sites)
#         # print(site_bond_scores)
#         # print(stoich_weighted_score, total_stoich)
#         # print("calculated score: {}".format(stoich_weighted_score / total_stoich))
#         return stoich_weighted_score / total_stoich
#
#
#     def get_site_bond_scores(self, bonds, shelx_file, n_bonds=4):
#         """
#             Determines priority for adding in site mixing.  Finds most problematic sites based on whether bonds
#             are shorter than expected.
#             TODO: What if bonds are longer than expected?
#             :param bonds: List of bond tuples (atom1, atom2, distance) for which to calculate priority
#             :return: List of (site number, site score) tuples sorted with decreasing priority.
#             """
#         bond_by_atom = defaultdict(lambda: [])
#         for bond in bonds:
#             bond_score = self.get_bond_score(bond, shelx_file)
#             bond_by_atom[bond[0]].append(bond_score)
#             bond_by_atom[bond[1]].append(bond_score)
#         # Average over scores from 4 shortest bonds
#         res = map(lambda tup: (tup[0], sum(sorted(tup[1])[:n_bonds])), bond_by_atom.items())
#         # Sort by which bonds are the shortest compared to what we'd expect
#         return sorted(res, key=lambda tup: tup[1])
#
#     def site_mixing_priority(self, bonds, shelx_file, n_bonds=2):
#
#         # Get associated site index
#         return map(lambda tup: (int(re.search("\d+", tup[0]).group(0)), tup[1]), self.get_site_bond_scores(bonds, shelx_file, n_bonds=n_bonds))
#
#     def get_bond_key(self, el1, el2):
#         return "-".join(sorted([el1, el2]))
#
#     def get_bond_length(self, el1, el2, shelx_file):
#         bond_key = self.get_bond_key(el1, el2)
#         if bond_key in self.bond_lengths:
#             return self.bond_lengths[bond_key]
#
#         if self.ml_model is not None:
#             return self.get_ml_prediction(el1, el2, shelx_file)
#
#         # Use pymatgen to get approximate bond length = sum of atomic radii
#         return Element(el1).atomic_radius + Element(el2).atomic_radius
#
#     def specie_to_el(self, specie_name):
#         return re.sub('\d', "", specie_name)
#
#     def get_ideal_bond_length(self, specie_name1, specie_name2, shelx_file):
#         """
#         Scores the likelihood of a bond based on its deviation from an ideal bond length.
#         :param specie_name1: One element in the bond
#         :param specie_name2: The other element in the bond
#         TODO: The ML model needs access to the .ins file to get the chemical system
#         If not, the ideal bond length is calculated as the sum of atomic radii.
#         TODO: This should be covalent radii instead
#         :return: ideal bond length
#         """
#         el1 = self.specie_to_el(specie_name1)
#         el2 = self.specie_to_el(specie_name2)
#         return self.get_bond_length(el1, el2, shelx_file)
#
#     def get_bonds(self, driver, shelx_file):
#         """
#         Get the bonds defined in the given shelx file
#
#         :param ins_file: SHELXFile object
#         :return res: List of bond tuples (element 1, element 2, bond length)
#         """
#         ins_file = copy.deepcopy(shelx_file)
#         ins_file.add_command("ACTA")
#         ins_file.remove_command("L.S.")
#         ins_file.add_command("L.S.", ["1"])
#         res = driver.run_SHELXTL(ins_file)
#         if res is None:
#             return []
#
#         with open(driver.cif_file) as f:
#             file_txt = f.read()
#             cif_file = CifParser.from_string(file_txt)
#
#         cif_dict = cif_file.as_dict().values()[0]
#         return zip(cif_dict["_geom_bond_atom_site_label_1"], cif_dict["_geom_bond_atom_site_label_2"],
#                    [float(x.replace("(", "").replace(")", "")) for x in cif_dict["_geom_bond_distance"]])
#
#     def get_n_missing_elements(self, shelx_file):
#         """
#         Get the number of missing elements based on the nominal formula
#         :param shelx_file: SHELX file
#         :return:
#         """
#         nominal_elements = set(map(lambda el: el.get_name(), shelx_file.elements))
#         result_elements = map(lambda site: site.el_string, shelx_file.get_all_sites())
#         return len(nominal_elements.difference(result_elements))
#
#     def get_stoichiometry_score(self, shelx_file):
#         """
#         Calculate a score for the given shelx file based on how well the composition agrees with the nominal formula
#         :param shelx_file: SHELX file
#         :return:
#         """
#         nominal_elements = set(map(lambda el: el.get_name(capitalize=True), shelx_file.elements))
#         nominal_formula = Composition(shelx_file.get_nominal_formula())
#         try:
#             analytic_formula = Composition(shelx_file.get_analytic_formula())
#         except:
#             print(shelx_file.get_analytic_formula())
#             quit()
#         score = 0.0
#         for el in nominal_elements:
#             diff = abs(nominal_formula.get_atomic_fraction(el) - analytic_formula.get_atomic_fraction(el))
#             if diff > 0.05:
#                 score += diff
#         return score
