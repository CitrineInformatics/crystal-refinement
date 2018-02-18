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
        if self.optimizer.log_output:
            print("Starting with {} sites".format(initial.res_file.get_n_sites()))
        if initial.r1 > self.optimizer.r1_threshold:
            self.reset_origin(initial)
            # self.try_add_q(initial)dy_bond_length(leaf)
        else:
            self.try_add_q(initial)
            for leaf in initial.get_leaves():
                self.try_remove_site(leaf)
        if self.optimizer.log_output:
            print("Finished with {} sites".format(", ".join([str(leaf.res_file.get_n_sites()) for leaf in initial.get_leaves()])))


    # def reset_origin(self, initial):
    #     new_ins = initial.get_res_copy()
    #     for i in reversed(range(1, len(new_ins.get_all_sites()))):
    #         new_ins.move_crystal_to_q(i)


    def identify_sites_by_bond_length(self, initial):
        """
        If the initial r1 is large, use bond lengths instead of the r1 score as the criteria for identifying which
        electron density peaks are atom sites
        """
        ins_file = initial.res_file
        shortest_possible_bond = self.optimizer.utils.get_shortest_bond(ins_file)
        prev_iteration = initial
        for i in range(5):
            new_ins = prev_iteration.get_res_copy()
            sites_to_remove = set()
            for bond in sorted(self.optimizer.utils.get_bonds(self.optimizer.driver, new_ins), key=lambda tup: tup[2]):
                # Threshold on how short the bonds are
                if bond[2] / shortest_possible_bond < 0.75:
                    a1_num = int(re.search('\d+', bond[0]).group(0))
                    a2_num = int(re.search('\d+', bond[1]).group(0))
                    index_to_remove = max(a1_num, a2_num)
                    for site in new_ins.get_sites_by_index(index_to_remove):
                        sites_to_remove.add(site)
            # no sites removed
            if len(sites_to_remove) == 0:
                break

            for site in sites_to_remove:
                new_ins.remove_site(site)

            annotation = "Deleted {} sites based on bond length".format(len(sites_to_remove))
            if self.optimizer.log_output:
                print("Trying step: {}".format(annotation))
            prev_iteration = self.optimizer.history.run_and_save(new_ins, prev_iteration, annotation)
            if prev_iteration is None:
                break


    def try_add_q(self, initial):
        """
        Try adding q peaks to main crystal sites if it decreases R value
        :return:
        """
        annotation = "Added q peak"
        base_ins = initial.get_res_copy()
        # This threshold could be scaled based on the potential atoms
        lightest_element = min([el.get_pymatgen_element().number for el in base_ins.elements])
        if base_ins.q_peaks[0].electron_density > lightest_element:
            base_ins.move_q_to_crystal()
            # Find best element for new site
            iterations = []
            for el in base_ins.elements:
                new_ins = base_ins.copy()
                new_ins.get_sites_by_index(new_ins.get_n_sites()).switch_element(el)
                if self.optimizer.log_output:
                    print("Trying step: {}".format(annotation))
                iteration = self.optimizer.history.run_iter(new_ins, initial, annotation)
                if iteration is not None:
                    iterations.append(iteration)
            iterations.sort(key=lambda i: i.r1)
            best_iter = iterations[0]
            # displacements = [x.displacement for x in best_iter.res_file.get_all_sites()]
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

        ins_file = initial.res_file
        r_penalty = 1.1
        bonds = self.optimizer.utils.get_bonds(self.optimizer.driver, ins_file)
        threshold = 0.1
        # print(len(ins_file.crystal_sites))
        while True:
            new_ins = initial.get_res_copy()
            sites_to_remove = set()
            for a1, a2, distance in bonds:
                ideal_distance = self.optimizer.utils.get_ideal_bond_length(a1, a2, new_ins)
                # if the distance is too small, remove the lower density atom
                if (ideal_distance - distance) / ideal_distance > threshold:
                    a1_num = int(re.search('\d+', a1).group(0))
                    a2_num = int(re.search('\d+', a2).group(0))
                    index_to_remove = max(a1_num, a2_num)
                    for site in new_ins.get_sites_by_index(index_to_remove):
                        sites_to_remove.add(site)
            # no sites removed
            if len(sites_to_remove) == 0:
                break

            for site in sites_to_remove:
                new_ins.remove_site(site)
            # ins_file.remove_sites_by_number(to_delete)
            # ins_file.renumber_sites()
            annotation = "Removed {} site(s)".format(len(sites_to_remove))
            if self.optimizer.log_output:
                print("Trying step: {}".format(annotation))
            cur_iter = self.optimizer.history.run_iter(new_ins, initial, annotation)

            if cur_iter is not None and cur_iter.r1 < initial.r1 * r_penalty:
                self.optimizer.history.save(cur_iter)
                break
            threshold *= 1.1
            # quit()


    def switch_elements(self, initial):
        ins_file = initial.res_file

        # Want to make changes from largest displacement to smallest
        sorted_sites = sorted(ins_file.get_all_sites(), key=lambda s: -s.displacement)

        for original_site in sorted_sites:
            for prev_iter in initial.get_leaves():
                iterations = []
                prev = original_site.get_name()
                for element in ins_file.elements:
                    new_ins = prev_iter.get_res_copy()
                    new_site = new_ins.get_sites_by_index(original_site.site_number)[0]
                    new_site.switch_element(element)
                    cur = new_site.get_name()
                    annotation = "Changed {} to {}".format(prev, cur)
                    if self.optimizer.log_output:
                        print("Trying step: {}".format(annotation))
                    iteration = self.optimizer.history.run_iter(new_ins, prev_iter, annotation)

                    if iteration is not None:
                        iterations.append(iteration)
                iterations.sort(key=lambda i: i.r1)
                # print(prev)
                for iteration in iterations:
                    # print(iteration.ins_file.crystal_sites[i].get_name(), [site.get_name() for site in iteration.ins_file.crystal_sites], iteration.r1, iteration.bond_score, iteration.get_score())
                    # if iteration.r1 - iterations[0].r1 < self.optimizer.r1_similarity_threshold:
                    # print(len(iterations))
                    # print([x.get_score() for x in iterations])
                    if iteration.r1 - iterations[0].r1 < self.optimizer.r1_similarity_threshold or \
                            iteration.get_score() / iterations[0].get_score() < self.optimizer.overall_score_ratio_threshold:
                        # print("saved, score ratio: {}".format(iteration.get_score() / iterations[0].get_score()))
                        self.optimizer.history.save(iteration)
            if original_site == sorted_sites[-1]:
                self.optimizer.history.clean_history(branch=initial)
            else:
                self.optimizer.history.clean_history(branch=initial, criteria="r1_only")
            # print("#"*50)
        # quit()


    def change_occupancy(self, initial):
        ins_file = initial.res_file

        # Want to make changes from largest displacement to smallest, but only in unmixed sites
        sorted_sites = sorted(ins_file.get_unmixed_sites(), key=lambda site: -site.displacement)
        displacements_list = map(lambda site: site.displacement, ins_file.get_all_sites())
        for original_site in sorted_sites:
            # don't change occupancy of mixed sites
            # if ins_file.crystal_sites[i].site_number in ins_file.mixed_site_numbers:
            #     continue
            for prev_iter in initial.get_leaves():
                new_ins = prev_iter.get_res_copy()
                # only change occupancy if displacement is >= 2 std deviations away from mean
                site_displacement = original_site.displacement
                filtered_displacements = filter(lambda d: d != site_displacement, displacements_list) # remove displacement of current site
                mean = np.mean(filtered_displacements)
                std = np.std(filtered_displacements)
                if abs(site_displacement - mean) / std < 2.0:
                    break

                for site in new_ins.get_sites_by_index(original_site.site_number):
                    new_ins.add_variable_occupancy(site)
                iteration = self.optimizer.history.run_iter(ins_file, initial, "Added variable occupancy for {}".format(original_site))
                new_site = iteration.res_file.get_sites_by_index(original_site.site_number)[0]
                # If changing the occupancy decreased r1, decreased the displacement, and resulted in an occupancy
                # that satisfies the threshold, add it to the history
                if iteration is not None \
                        and iteration.r1 < prev_iter.r1 \
                        and float(iteration.res_file.fvar_vals[-1]) < (1 - self.optimizer.occupancy_threshold) \
                        and new_site.displacement < site_displacement:
                    self.optimizer.history.save(iteration)


    def try_site_mixing(self, initial):
        """
        Try allowing site mixing, where a single crystal site might have mixed occupancy between two different elements
        Only handles pair-wise mixing.
        :param driver: SHELX driver
        """

        # pairs = self.optimizer.utils.mixing_pairs
        pairs = self.optimizer.utils.mixing_pairs

        # Keep adding site mixing until we've already tried adding site mixing at the top priority sites
        tried = set()
        self.do_site_mixing(initial, tried, pairs)


    def do_site_mixing(self, initial, tried, pairs):
        ins_file = initial.res_file
        bonds = self.optimizer.utils.get_bonds(self.optimizer.driver, ins_file)
        mixing_priority = self.optimizer.utils.site_mixing_priority(bonds, ins_file)
        # In case of ties, find all top tied priorities
        top_priority_score = mixing_priority[0][1]
        top_priority = [priority[0] for priority in mixing_priority if priority[1] - top_priority_score < 0.5]
        # For each of these top priorities, try site mixing
        if all([i in tried for i in top_priority]):
            if self.optimizer.ensure_identified_elements:
                self.ensure_identified_elements_are_present(initial)
            return
        for i in top_priority:
            if i in tried:
                continue
            for prev_iter in initial.get_leaves():
                tried.add(i)
                iterations = []
                for pair in pairs:
                    new_ins = prev_iter.get_res_copy()
                    new_ins.add_equal_site_mixing(i, pair)
                    annotation = "Mixing {} and {} on site {} equally".format(pair[0].get_name(), pair[1].get_name(), i)
                    if self.optimizer.log_output:
                        print("Trying step: {}".format(annotation))
                    iteration = self.optimizer.history.run_iter(new_ins, prev_iter, annotation)

                    if iteration is not None:
                        if prev_iter.r1 - iteration.r1 < self.optimizer.r1_similarity_threshold:
                            new_ins = iteration.get_res_copy()
                            new_ins.add_site_mixing_variable_occupancy(i)
                            annotation = "Adding variable occupancy for {} and {} on site {}".format(pair[0].get_name(), pair[1].get_name(), i)
                            if self.optimizer.log_output:
                                print("Trying step: {}".format(annotation))
                            variable_occupancy_iteration = self.optimizer.history.run_iter(new_ins, iteration, annotation)
                            if variable_occupancy_iteration is not None:
                                occupancy_var = float(variable_occupancy_iteration.get_res_copy().fvar_vals[-1])
                                # Only include occupancies that are actually split
                                if occupancy_var > self.optimizer.occupancy_threshold and occupancy_var < (1 - self.optimizer.occupancy_threshold):
                                    iterations.append(variable_occupancy_iteration)

                if len(iterations) == 0:
                    continue
                iterations.sort(key=lambda i: i.r1)
                best_r1 = min([prev_iter.r1, iterations[0].r1])
                if prev_iter.r1 - best_r1 < self.optimizer.r1_similarity_threshold:
                    prev_iter.propagate()
                if iterations[0].r1 - best_r1 < self.optimizer.r1_similarity_threshold:
                    self.optimizer.history.save(iterations[0])
        self.optimizer.history.clean_history(branch=initial)
        if self.optimizer.log_output:
            print("Cleared branch history for {} site path".format(len(ins_file.get_all_sites())))
        for leaf in initial.get_leaves():
            self.do_site_mixing(leaf, tried.union(set(top_priority)), pairs)


    def try_anisotropy(self, initial):
        """
        Test if adding anisotropy reduces R1 value.  If it does, do so.
        :return:
        """

        new_ins = initial.get_res_copy()

        #  Try with anisotropy
        new_ins.add_anisotropy()
        annotation = "Added anisotropy"
        if self.optimizer.log_output:
            print("Trying step: {}".format(annotation))
        iteration = self.optimizer.history.run_iter(new_ins, initial, annotation)

        if iteration is not None:
            self.optimizer.history.save(iteration)
            if iteration.r1 > initial.r1:
                initial.propagate()


    def try_exti(self, initial):
        """
        Test if adding extinguishing reduces R1 value.  If it does, do so.
        :return:
        """
        new_ins = initial.get_res_copy()

        #  Try with extinguishing
        new_ins.add_exti()
        annotation = "Added extinction"
        if self.optimizer.log_output:
            print("Trying step: {}".format(annotation))
        iteration = self.optimizer.history.run_iter(new_ins, initial, "Added extinction")


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
        new_ins = initial.get_res_copy()

        #  Try with suggested
        new_ins.remove_command("WGHT")
        new_ins.commands.append(("WGHT", new_ins.suggested_weight_vals))
        annotation = "Used suggested weights"
        if self.optimizer.log_output:
            print("Trying step: {}".format(annotation))
        iteration = self.optimizer.history.run_iter(new_ins, initial, "Used suggested weights")

        if iteration is not None:
            self.optimizer.history.save(iteration)
            if iteration.r1 > initial.r1:
                initial.propagate()


    def ensure_identified_elements_are_present(self, initial):
        """
        Filter the given branch if there are missing elements
        :param initial:
        :return:
        """

        new_ins = initial.get_res_copy()

        identified_elements = set(map(lambda el: el.get_name(), new_ins.elements))
        result_elements = set()
        for site in new_ins.get_all_sites():
            result_elements.add(site.el_string.upper())
        missing = identified_elements.difference(result_elements)
        if len(missing) > 0:
            if self.optimizer.log_output:
                print("Truncated branch due to missing elements: {}".format(" ".join(missing)))
            self.optimizer.history.clean_history(n_to_keep=0, branch=initial)