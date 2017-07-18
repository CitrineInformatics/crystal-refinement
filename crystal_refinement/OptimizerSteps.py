import itertools
import re

import numpy as np
from pymatgen.core.composition import Element



class OptimizerSteps:
    """
    Class for performing single crystal refinement
    """

    def __init__(self, optimizer):
        self.optimizer = optimizer


    def identify_sites(self, initial):
        if initial.r1 > self.optimizer.r1_threshold:
            self.try_add_q(initial)
            self.identify_sites_by_bond_length(initial)
        else:
            self.try_add_q(initial)
            for leaf in initial.get_leaves():
                self.try_remove_site(leaf)


    def identify_sites_by_bond_length(self, initial):
        """
        If the initial r1 is large, use bond lengths instead of the r1 score as the criteria for identifying which
        electron density peaks are atom sites
        """
        ins_file = initial.get_res()
        shortest_possible_bond = self.optimizer.utils.get_shortest_bond(ins_file)
        prev_iteration = initial
        for i in range(5):
            ins_file = prev_iteration.get_res()
            to_delete = set()
            for bond in sorted(self.optimizer.utils.get_bonds(self.optimizer.driver, ins_file), key=lambda tup: tup[2]):
                # Threshold on how short the bonds are
                if bond[2] / shortest_possible_bond < 0.75:
                    a1_num = int(re.search('\d+', bond[0]).group(0))
                    a2_num = int(re.search('\d+', bond[1]).group(0))
                    to_delete.add(max(a1_num, a2_num))
            if len(to_delete) == 0:
                break
            ins_file.remove_sites_by_number(to_delete)
            ins_file.renumber_sites()
            prev_iteration = self.optimizer.history.run_and_save(ins_file, prev_iteration)


    def try_add_q(self, initial):
        """
        Try adding q peaks to main crystal sites if it decreases R value
        :return:
        """
        annotation = "Added q peak"
        ins_file = initial.get_res()
        # This threshold could be scaled based on the potential atoms
        lightest_element = min([Element(el.capitalize()).number for el in ins_file.elements])
        if ins_file.q_peaks[0].electron_density > lightest_element:
            ins_file.move_q_to_crystal()
            # Find best element for new site
            num_elems = len(ins_file.elements)
            iterations = []
            for elem in range(1, num_elems + 1):
                ins_file.change_element(len(ins_file.crystal_sites) - 1, elem)
                iteration = self.optimizer.history.run_iter(ins_file, initial, annotation)
                if iteration is not None:
                    iterations.append(iteration)
            iterations.sort(key=lambda i: i.r1)
            best_iter = iterations[0]
            displacements = [x.displacement for x in best_iter.res_file.crystal_sites]
            if best_iter.r1 < initial.r1:  # and displacements[-1] < (np.mean(displacements[:-1]) + 2.0 * np.std(displacements[:-1])):
                self.optimizer.history.save([iterations[0]])
                for leaf in initial.get_leaves():
                    #   If adding one peak helped, recursively try adding another peak until it stops helping
                    self.try_add_q(leaf)


    def try_remove_site(self, initial):
        """
        Remove crystal sites if they result in bond distances that are too short
        :param ml_model: if True, a machine learning model is used to predict the correct bond length.
                         if False, the ideal bond length is approximated as the sum of the atomic radii
        :return:
        """

        ins_file = initial.get_res()
        r_penalty = 1.1
        bonds = self.optimizer.utils.get_bonds(self.optimizer.driver, ins_file)
        threshold = 0.1
        # print(len(ins_file.crystal_sites))
        while True:
            ins_file = initial.get_res()
            to_delete = set()
            for a1, a2, distance in bonds:
                ideal_distance = self.optimizer.utils.get_ideal_bond_length(a1, a2, ins_file)
                # if the distance is too small, remove the lower density atom
                if (ideal_distance - distance) / ideal_distance > threshold:
                    a1_num = int(re.search('\d+', a1).group(0))
                    a2_num = int(re.search('\d+', a2).group(0))
                    to_delete.add(max(a1_num, a2_num))
            # no sites removed
            if len(to_delete) == 0:
                break
            ins_file.remove_sites_by_number(to_delete)
            ins_file.renumber_sites()

            cur_iter = self.optimizer.history.run_iter(ins_file, initial, "Removed {} site(s)".format(len(to_delete)))

            if cur_iter is not None and cur_iter.r1 < initial.r1 * r_penalty:
                self.optimizer.history.save(cur_iter)
                break
            threshold *= 1.1
            # quit()


    def switch_elements(self, initial):
        ins_file = initial.get_res()

        # Want to make changes from largest displacement to smallest
        sorted_sites = sorted(ins_file.crystal_sites, key=lambda s: -s.displacement)
        order = [s.site_number - 1 for s in sorted_sites]
        num_elems = len(ins_file.elements)

        # print(ins_file.get_crystal_sites_text())
        for i in order:
            for prev_iter in initial.get_leaves():
                ins_file = prev_iter.get_res()
                iterations = []
                for elem in range(1, num_elems + 1):
                    ins_file.change_element(i, elem)
                    prev = prev_iter.res_file.crystal_sites[i].get_name()
                    cur = ins_file.crystal_sites[i].get_name()
                    iteration = self.optimizer.history.run_iter(ins_file, prev_iter, "Changed {} to {}".format(prev, cur))

                    if iteration is not None:
                        iterations.append(iteration)
                iterations.sort(key=lambda i: i.r1)
                # print(prev_iter.res_file.crystal_sites[i].get_name())
                for iteration in iterations:
                    # print(iteration.ins_file.crystal_sites[i].get_name(), iteration.r1, iteration.bond_score, iteration.get_score())
                    if iteration.r1 - iterations[0].r1 < self.optimizer.r1_similarity_threshold:
                        # print("saved")
                        self.optimizer.history.save(iteration)
            self.optimizer.history.clean_history()
            # print("#"*50)
        # quit()


    def change_occupancy(self, initial):
        ins_file = initial.get_res()

        # Want to make changes from largest displacement to smallest
        displacements = list(map((lambda x: x.displacement), ins_file.crystal_sites))

        for i, displacement in sorted(enumerate(displacements), key=lambda tup: -tup[1]):
            # don't change occupancy of mixed sites
            if ins_file.crystal_sites[i].site_number in ins_file.mixed_site_numbers:
                continue
            for prev_iter in initial.get_leaves():
                ins_file = prev_iter.get_res()

                # only change occupancy if displacement is >= 2 std deviations away from mean
                mean = np.mean(displacements[:i] + displacements[i + 1:])
                std = np.std(displacements[:i] + displacements[i + 1:])
                if abs(displacement - mean) / std < 2.0:
                    break

                ins_file.add_variable_occupancy(i)

                site = ins_file.crystal_sites[i].get_name()
                iteration = self.optimizer.history.run_iter(ins_file, initial, "Added variable occupancy for {}".format(site))


                # If changing the occupancy decreased r1, decreased the displacement, and resulted in an occupancy
                # that satisfies the threshold, add it to the history
                if iteration is not None \
                        and iteration.r1 < prev_iter.r1 \
                        and float(iteration.res_file.fvar_vals[-1]) < (1 - self.optimizer.occupancy_threshold) \
                        and iteration.res_file.crystal_sites[i].displacement < displacement:
                    self.optimizer.history.save(iteration)


    def try_site_mixing(self, initial):
        """
        Try allowing site mixing, where a single crystal site might have mixed occupancy between two different elements
        Only handles pair-wise mixing.
        :param driver: SHELX driver
        """

        pairs = self.optimizer.utils.mixing_pairs

        # Keep adding site mixing until we've already tried adding site mixing at the top priority sites
        tried = set()
        self.do_site_mixing(initial, tried, pairs)


    def do_site_mixing(self, initial, tried, pairs):
        ins_file = initial.get_res()
        bonds = self.optimizer.utils.get_bonds(self.optimizer.driver, ins_file)
        mixing_priority = list(self.optimizer.utils.site_mixing_priority(bonds, ins_file))
        # In case of ties, find all top tied priorities
        top_priority_score = mixing_priority[0][1]
        top_priority = [priority[0] for priority in mixing_priority if priority[1] - top_priority_score < 0.5]
        # For each of these top priorities, try site mixing
        if all([i in tried for i in top_priority]):
            return
        for i in top_priority:
            if i in tried:
                continue
            for prev_iter in initial.get_leaves():
                tried.add(i)
                iterations = []
                for pair in pairs:
                    ins_file = prev_iter.get_res()
                    ins_file.add_site_mixing(site_number=i, mixing_element_indices=pair)

                    # Graph annotation
                    mix = "{} and {}".format(ins_file.elements[pair[0]], ins_file.elements[pair[1]])
                    iteration = self.optimizer.history.run_iter(ins_file, prev_iter, "Mixing {} on site {}".format(mix, i))

                    if iteration is not None:
                        occupancy_var = float(iteration.res_file.fvar_vals[-1])
                        # Only include occupancies that are actually split
                        if occupancy_var > self.optimizer.occupancy_threshold and occupancy_var < (1 - self.optimizer.occupancy_threshold):
                            iterations.append(iteration)

                if len(iterations) == 0:
                    continue
                iterations.sort(key=lambda i: i.r1)
                best_r1 = min([prev_iter.r1, iterations[0].r1])
                if prev_iter.r1 - best_r1 < self.optimizer.r1_similarity_threshold:
                    prev_iter.propagate()
                if iterations[0].r1 - best_r1 < self.optimizer.r1_similarity_threshold:
                    self.optimizer.history.save(iterations[0])
        self.optimizer.history.clean_history()
        for leaf in initial.get_leaves():
            self.do_site_mixing(leaf, tried.union(set(top_priority)), pairs)


    def try_anisotropy(self, initial):
        """
        Test if adding anisotropy reduces R1 value.  If it does, do so.
        :return:
        """

        ins_file = initial.get_res()

        #  Try with anisotropy
        ins_file.add_anisotropy()

        iteration = self.optimizer.history.run_iter(ins_file, initial, "Added anisotropy")

        if iteration is not None:
            self.optimizer.history.save(iteration)
            if iteration.r1 > initial.r1:
                initial.propagate()


    def try_exti(self, initial):
        """
        Test if adding extinguishing reduces R1 value.  If it does, do so.
        :return:
        """
        ins_file = initial.get_res()

        #  Try with extinguishing
        ins_file.add_exti()

        iteration = self.optimizer.history.run_iter(ins_file, initial, "Added extinction")


        # If exti helped, add it to the history
        if iteration is not None:
            self.optimizer.history.save(iteration)
            if iteration.r1 > initial.r1:
                initial.propagate()


    def use_suggested_weights(self, initial):
        """
        Stop re-initializing weights each time--use previously suggested weights
        :return:
        """
        ins_file = initial.get_res()

        #  Try with suggested
        ins_file.remove_command("WGHT")
        ins_file.commands.append(("WGHT", ins_file.suggested_weight_vals))

        iteration = self.optimizer.history.run_iter(ins_file, initial, "Used suggested weights")

        if iteration is not None:
            self.optimizer.history.save(iteration)
            if iteration.r1 > initial.r1:
                initial.propagate()