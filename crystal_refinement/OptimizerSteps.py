def identify_sites(self, initial):
    if initial.r1 > self.r1_threshold:
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
    shortest_possible_bond = self.get_shortest_bond(ins_file)
    prev_iteration = initial
    for i in range(5):
        ins_file = prev_iteration.get_res()
        to_delete = set()
        for bond in sorted(self.get_bonds(ins_file), key=lambda tup: tup[2]):
            # Threshold on how short the bonds are
            if bond[2] / shortest_possible_bond < 0.5:
                a1_num = int(re.search('\d+', bond[0]).group(0))
                a2_num = int(re.search('\d+', bond[1]).group(0))
                to_delete.add(max(a1_num, a2_num))
        if len(to_delete) == 0:
            break
        ins_file.remove_sites_by_number(to_delete)
        prev_iteration = self.history.run_and_save(ins_file, prev_iteration)

def try_add_q(self, initial):
    """
    Try adding q peaks to main crystal sites if it decreases R value
    :return:
    """

    ins_file = initial.get_res()

    # This threshold could be scaled based on the potential atoms
    if ins_file.q_peaks[0].electron_density > 50:
        ins_file.move_q_to_crystal()
        # Find best element for new site
        num_elems = len(ins_file.elements)
        iterations = []
        for elem in range(1, num_elems + 1):
            ins_file.change_element(len(ins_file.crystal_sites)-1, elem)
            iterations.append(self.history.run_iter(ins_file, initial))
        iterations.sort(key=lambda i: i.r1)
        if iterations[0].r1 - initial.r1 < self.r1_similarity_threshold:
            self.history.save([iterations[0]])
            for leaf in initial.get_leaves():
                #   If adding one peak helped, recursively try adding another peak until it stops helping
                self.try_add_q(leaf)

def try_remove_site(self, initial, use_ml_model=False):
    """
    Remove crystal sites if they result in bond distances that are too short
    :param ml_model: if True, a machine learning model is used to predict the correct bond length.
                     if False, the ideal bond length is approximated as the sum of the atomic radii
    :return:
    """

    ins_file = initial.get_res()

    r_penalty = 1.1
    bonds = self.get_bonds(ins_file)
    threshold = 0.1
    while True:
        ins_file = initial.get_res()
        to_delete = set()
        for a1, a2, distance in bonds:
            ideal_distance = utils.get_ideal_bond_length(a1, a2, use_ml_model)
            # if the distance is too small, remove the lower density atom
            if (ideal_distance - distance) / ideal_distance > threshold:
                a1_num = int(re.search('\d+', a1).group(0))
                a2_num = int(re.search('\d+', a2).group(0))
                to_delete.add(max(a1_num, a2_num))
        # no sites removed
        if len(to_delete) == 0:
            break
        ins_file.remove_sites_by_number(to_delete)
        cur_iter = self.history.run_iter(ins_file, initial)
        if cur_iter.r1 < initial.r1 * r_penalty:
            self.history.save(cur_iter)
            break
        threshold *= 1.1


    def try_site_mixing(self, initial):
        """
        Try allowing site mixing, where a single crystal site might have mixed occupancy between two different elements
        Only handles pair-wise mixing.
        :param driver: SHELX driver
        """
        ins_file = initial.get_res()
        element_list = [Element(el.capitalize()) for el in ins_file.elements]
        pairs = []
        probability_threshold = 2E-4

        # For all elements in compound, calculate substitution probabilities
        # If substitution probability is > probability_threshold, then save it to pairs list.
        for i1, i2 in itertools.combinations(range(len(element_list)), 2):
            e1 = element_list[i1]
            e2 = element_list[i2]
            sp = utils.get_substitution_probability(e1, e2)
            if sp > probability_threshold:
                pairs.append(([i1, i2], sp))

        # Sort pairs by substitution probability (largest to smallest)
        pairs = [tup[0] for tup in sorted(pairs, key=lambda tup: -tup[1])]

        ins_file = initial.get_res()
        bonds = self.get_bonds(ins_file)
        mixing_priority = [x[0] for x in utils.site_mixing_priority(bonds)]
        for i in mixing_priority:
            for prev_iter in initial.get_leaves():
                iterations = []
                for pair in pairs:
                    ins_file = prev_iter.get_res()
                    ins_file.add_site_mixing(site_number=i, mixing_element_indices=pair)
                    iteration = self.history.run_iter(ins_file, prev_iter)
                    if iteration is not None:
                        occupancy_var = float(iteration.res_file.fvar_vals[-1])
                        # Only include occupancies that are actually split
                        if occupancy_var > 0.02 and occupancy_var < 0.98:
                            iterations.append(iteration)
                if len(iterations) == 0:
                    continue
                best_r1 = min([prev_iter.r1] + [x.r1 for x in iterations])
                if prev_iter.r1 - best_r1 < self.r1_similarity_threshold:
                    prev_iter.propagate()
                for iteration in iterations:
                    if iteration.r1 - best_r1 < self.r1_similarity_threshold:
                        self.history.save(iteration)


    # def try_site_mixing(self, initial):
    #     """
    #     Try allowing site mixing, where a single crystal site might have mixed occupancy between two different elements
    #     Only handles pair-wise mixing.
    #     :param driver: SHELX driver
    #     """
    #     ins_file = initial.get_res()
    #     element_list = [Element(el.capitalize()) for el in ins_file.elements]
    #     pairs = []
    #     probability_threshold = 2E-4
    #
    #     # For all elements in compound, calculate substitution probabilities
    #     # If substitution probability is > probability_threshold, then save it to pairs list.
    #     for i1, i2 in itertools.combinations(range(len(element_list)), 2):
    #         e1 = element_list[i1]
    #         e2 = element_list[i2]
    #         sp = utils.get_substitution_probability(e1, e2)
    #         if sp > probability_threshold:
    #             pairs.append(([i1, i2], sp))
    #
    #     # Sort pairs by substitution probability (largest to smallest)
    #     pairs = [tup[0] for tup in sorted(pairs, key=lambda tup: -tup[1])]
    #
    #     # Keep adding site mixing until we've already tried adding site mixing at the top priority sites
    #     tried = set()
    #     self.do_site_mixing(initial, tried, pairs)
    #
    # def do_site_mixing(self, initial, tried, pairs):
    #     ins_file = initial.get_res()
    #     bonds = self.get_bonds(ins_file)
    #     mixing_priority = utils.site_mixing_priority(bonds)
    #     # In case of ties, find all top tied priorities
    #     top_priority_score = mixing_priority[0][1]
    #     top_priority = [priority[0] for priority in mixing_priority if np.abs(priority[1] - top_priority_score) < 0.001]
    #
    #     # For each of these top priorities, try site mixing
    #     if all([i in tried for i in top_priority]):
    #         return
    #     for i in top_priority:
    #         if i in tried:
    #             continue
    #         for prev_iter in initial.get_leaves():
    #             tried.add(i)
    #             iterations = []
    #             for pair in pairs:
    #                 ins_file = prev_iter.get_res()
    #                 ins_file.add_site_mixing(site_number=i, mixing_element_indices=pair)
    #                 iteration = self.history.run_iter(ins_file, prev_iter)
    #                 if iteration is not None:
    #                     occupancy_var = float(iteration.res_file.fvar_vals[-1])
    #                     # Only include occupancies that are actually split
    #                     if occupancy_var > 0.02 and occupancy_var < 0.98:
    #                         iterations.append(iteration)
    #             if len(iterations) == 0:
    #                 continue
    #             iterations.sort(key=lambda i: i.r1)
    #             if iterations[0].r1 < prev_iter.r1:
    #                 self.history.save(iterations[0])
    #                 self.do_site_mixing(iterations[0], tried.union(set(top_priority)), pairs)