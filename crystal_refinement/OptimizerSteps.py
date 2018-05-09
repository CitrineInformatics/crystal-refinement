import re
import numpy as np
import crystal_refinement.utils.bond_utils as bond_utils


def identify_sites(initial, optimizer):
    """
    Identify which electron density peaks are atom sites
    """
    if optimizer.log_output:
        print("Starting with {} sites".format(initial.res_file.get_n_sites()))
    # Handle cases with large initial r1
    if initial.r1 > optimizer.r1_threshold:
        if initial.res_file.get_n_sites() == 1 and initial.res_file.is_small_cubic_structure():
            reset_origin(initial, optimizer)
        else:
            identify_sites_by_bond_length(initial, optimizer)
    else:
        try_add_q(initial, optimizer)
        for leaf in initial.get_leaves():
            try_remove_site(leaf, optimizer)
    if optimizer.log_output:
        print("Finished with {} sites".format(", ".join([str(leaf.res_file.get_n_sites()) for leaf in initial.get_leaves()])))


def reset_origin(initial, optimizer):
    """
    Performs an origin reset for the case where there is 1 site. Moves the single site to 0, 0, 0 and sets the
    site occupancy and displacement to appropriate values for a small cubic structure. Then tries to add q peaks
    as additional sites.
    """
    new_ins = initial.get_res_copy()
    for i in reversed(range(1, new_ins.get_n_sites())):
        new_ins.move_crystal_to_q(i)
    for site in new_ins.get_sites_by_index(1):
        site.position = np.asarray([0.0, 0.0, 0.0])
        site.occupancy = 0.02083
        site.displacement = 0.05
    iteration = optimizer.history.run_iter(new_ins, initial, "Reset origin")
    if iteration is not None:
        optimizer.history.save(iteration)
        # Try to add additional sites
        try_add_q(iteration, optimizer)


def identify_sites_by_bond_length(initial, optimizer):
    """
    If the initial r1 is large, use bond lengths instead of the r1 score as the criteria for identifying which
    electron density peaks are atom sites
    """

    # Identify shortest possible bond based on available elements
    shortest_possible_bond = optimizer.cache.get_shortest_bond()
    prev_iteration = initial
    for i in range(5):
        new_ins = prev_iteration.get_res_copy()
        sites_to_remove = set()
        for bond in sorted(bond_utils.get_bonds(optimizer.driver, new_ins)):
            # Threshold on how short the bonds are
            if bond.length / shortest_possible_bond < 0.75:
                a1_num = int(re.search('\d+', bond.el1).group(0))
                a2_num = int(re.search('\d+', bond.el2).group(0))
                index_to_remove = max(a1_num, a2_num)
                for site in new_ins.get_sites_by_index(index_to_remove):
                    sites_to_remove.add(site)
        # no sites removed
        if len(sites_to_remove) == 0:
            if i == 0:
                try_add_q(initial, optimizer)
            break

        for site in sites_to_remove:
            new_ins.remove_site(site)

        annotation = "Deleted {} sites based on bond length".format(len(sites_to_remove))
        if optimizer.log_output:
            print("Trying step: {}".format(annotation))
        prev_iteration = optimizer.history.run_and_save(new_ins, prev_iteration, annotation)
        if prev_iteration is None:
            break


def try_add_q(initial, optimizer):
    """
    Try adding q peaks to main crystal sites if it decreases R value
    """
    base_ins = initial.get_res_copy()
    lightest_element = min([el.get_pymatgen_element().number for el in base_ins.elements])
    if base_ins.q_peaks[0].electron_density > lightest_element:
        annotation = "Added q peak {}".format(base_ins.q_peaks[0].to_string())
        base_ins.move_q_to_crystal()
        # Find best element for new site
        iterations = []
        for el in base_ins.elements:
            new_ins = base_ins.copy()
            for site in new_ins.get_sites_by_index(new_ins.get_n_sites()):
                site.switch_element(el)
            if optimizer.log_output:
                print("Trying step: {}".format(annotation))
            iteration = optimizer.history.run_iter(new_ins, initial, annotation)
            if iteration is not None:
                iterations.append(iteration)
        iterations.sort(key=lambda i: i.r1)
        best_iter = iterations[0]
        if best_iter.r1 < initial.r1:
            optimizer.history.save([iterations[0]])
            for leaf in initial.get_leaves():
                #   If adding one peak helped, recursively try adding another peak until it stops helping
                try_add_q(leaf, optimizer)


def try_remove_site(initial, optimizer):
    """
    Remove crystal sites if they result in bond distances that are too short
    """
    ins_file = initial.res_file
    r_penalty = 1.1
    bonds = bond_utils.get_bonds(optimizer.driver, ins_file)
    threshold = 0.05
    removal_sets_tried = set()
    while True:
        new_ins = initial.get_res_copy()
        sites_to_remove = set()
        distance_ratios = []
        for bond in bonds:
            ideal_distance = optimizer.cache.get_ideal_bond_length(*bond.get_normalized_element_names())
            distance_ratio = (ideal_distance - bond.length) / ideal_distance
            a1_num = int(re.search('\d+', bond.el1).group(0))
            a2_num = int(re.search('\d+', bond.el2).group(0))
            index_to_remove = max(a1_num, a2_num)
            for site in new_ins.get_sites_by_index(index_to_remove):
                distance_ratios.append((distance_ratio, site))
            # if the distance is too small, remove the site with lower electron density
            if (ideal_distance - bond.length) / ideal_distance > threshold:
                for site in new_ins.get_sites_by_index(index_to_remove):
                    sites_to_remove.add(site)
        set_code = ", ".join(sorted(map(lambda site: site.get_name(), sites_to_remove)))
        if set_code in removal_sets_tried:
            threshold *= 1.1
            continue
        else:
            removal_sets_tried.add(set_code)
        # no sites removed
        if len(sites_to_remove) == 0:
            break
        new_ins.remove_sites(sites_to_remove)

        annotation = "Removed sites {}".format(set_code)
        if optimizer.log_output:
            print("Trying step: {}".format(annotation))
        cur_iter = optimizer.history.run_iter(new_ins, initial, annotation)
        if cur_iter is not None and cur_iter.r1 < initial.r1 * r_penalty:
            optimizer.history.save(cur_iter)
            break
        threshold *= 1.1


def switch_elements(initial, optimizer):
    """
    Switch elements to find the correct elemental site assignments.
    """
    ins_file = initial.res_file

    # Want to make changes from largest displacement to smallest
    sorted_sites = sorted(ins_file.get_all_sites(), key=lambda s: -s.displacement)

    for original_site in sorted_sites:
        # Try switching the selected site in every path generated so far
        for prev_iter in initial.get_leaves():
            prev = original_site.get_name()
            # Try switching each available element in
            for element in ins_file.elements:
                new_ins = prev_iter.get_res_copy()
                new_site = new_ins.get_sites_by_index(original_site.site_number)[0]
                new_site.switch_element(element)
                cur = new_site.get_name()
                annotation = "Changed {} to {}".format(prev, cur)
                if optimizer.log_output:
                    print("Trying step: {}".format(annotation))
                iteration = optimizer.history.run_iter(new_ins, prev_iter, annotation)

                # If the iteration didn't fail, save it
                if iteration is not None:
                    optimizer.history.save(iteration)

        # Pick the best paths after all switches have been made for the selected site
        if original_site == sorted_sites[-1]:
            # Criteria is overall score and r1 score in the last iteration
            optimizer.history.clean_history(branch=initial, criteria=["overall_score", "r1_only"])
        else:
            # Criteria is only r1 score in the last iteration
            optimizer.history.clean_history(branch=initial, criteria="r1_only")


def change_occupancy(initial, optimizer):
    """
    Try partial occupancy for the sites.
    """
    ins_file = initial.res_file
    # Want to make changes from largest displacement to smallest, but only in unmixed sites
    sorted_sites = sorted(ins_file.get_unmixed_sites(), key=lambda site: -site.displacement)
    displacements_list = map(lambda site: site.displacement, ins_file.get_all_sites())
    for original_site in sorted_sites:
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
            iteration = optimizer.history.run_iter(new_ins, initial, "Added variable occupancy for {}".format(original_site.get_name()))
            # If changing the occupancy decreased r1, decreased the displacement, and resulted in an occupancy
            # that satisfies the threshold, add it to the history
            if iteration is not None:
                new_site = iteration.res_file.get_sites_by_index(original_site.site_number)[0]
                if iteration.r1 < prev_iter.r1 \
                    and float(iteration.res_file.fvar_vals[-1]) < (1 - optimizer.occupancy_threshold) \
                    and new_site.displacement < site_displacement:
                    optimizer.history.save(iteration)


def try_site_mixing(initial, optimizer):
    """
    Try allowing site mixing, where a single crystal site might have mixed occupancy between two different elements.
    Only handles pair-wise mixing.
    """

    pairs = optimizer.cache.mixing_pairs

    # Set of sites where mixing has already been applied
    tried = set()
    _do_site_mixing(initial, tried, pairs, optimizer)


def _do_site_mixing(initial, tried, pairs, optimizer):
    """
    Actually try site mixing.
    This function identifies top priority sites for mixing. For each mixing pair, it tries equal ratio mixing then
    allows the ratio to vary freely. The function then recursively continues that process. It will stop when the
    top priority sites have all already been tried.
    :param initial: Starting .ins file
    :param tried: the sites where mixing has already been tried
    :param pairs: pairs of elements that are allowed to mix
    :return:
    """

    ins_file = initial.res_file
    bonds = bond_utils.get_bonds(optimizer.driver, ins_file)

    site_bond_scores = bond_utils.get_site_bond_scores(bonds, optimizer.cache, ins_file, n_bonds=2)
    mixing_priority = map(lambda tup: (int(re.search("\d+", tup[0]).group(0)), tup[1]), site_bond_scores)

    # In case of ties, find all top tied priorities
    top_priority_score = mixing_priority[0][1]
    top_priority = [priority[0] for priority in mixing_priority if priority[1] - top_priority_score < 0.5]

    # Check if we have tried mixing on all the top priority sites already. If so, stop recursion.
    if all([i in tried for i in top_priority]):
        if optimizer.ensure_identified_elements:
            ensure_identified_elements_are_present(initial, optimizer)
        return

    # For each of these top priorities, try site mixing
    for i in top_priority:
        if i in tried:
            continue
        for prev_iter in initial.get_leaves():
            tried.add(i)
            iterations = []

            # Try all the mixing pairs on the selected site
            for pair in pairs:
                new_ins = prev_iter.get_res_copy()

                # Start with equal mixing on the site
                new_ins.add_equal_site_mixing(i, pair)
                annotation = "Mixing {} and {} on site {} equally".format(pair[0].get_name(), pair[1].get_name(), i)
                if optimizer.log_output:
                    print("Trying step: {}".format(annotation))
                iteration = optimizer.history.run_iter(new_ins, prev_iter, annotation)

                # Check if the equal mixing was successful based on validity and r1 score
                if iteration is not None and iteration.r1 - prev_iter.r1 < optimizer.r1_similarity_threshold:
                    optimizer.history.save(iteration)
                    if prev_iter.r1 - iteration.r1 < optimizer.r1_similarity_threshold:
                        prev_iter.propagate()
                    new_ins = iteration.get_res_copy()

                    # Try allowing the mixing ratio to be variable
                    new_ins.add_site_mixing_variable_occupancy(i)
                    annotation = "Adding variable occupancy for {} and {} on site {}".format(pair[0].get_name(), pair[1].get_name(), i)
                    if optimizer.log_output:
                        print("Trying step: {}".format(annotation))
                    variable_occupancy_iteration = optimizer.history.run_iter(new_ins, iteration, annotation)

                    # Check if the variable mixing ratio was successful
                    if variable_occupancy_iteration is not None:
                        occupancy_var = float(variable_occupancy_iteration.get_res_copy().fvar_vals[-1])
                        if occupancy_var > optimizer.occupancy_threshold and occupancy_var < (1 - optimizer.occupancy_threshold):
                            iterations.append(variable_occupancy_iteration)

            # If no element pairs were successful, move on
            if len(iterations) == 0:
                continue

            # Pick the best element pair of all tried
            iterations.sort(key=lambda x: x.r1)
            best_r1 = min([prev_iter.r1, iterations[0].r1])
            if iterations[0].r1 - best_r1 < optimizer.r1_similarity_threshold:
                optimizer.history.save(iterations[0])

    # Prune the branches for this set of top priority mixing sites
    optimizer.history.clean_history(branch=initial, criteria=["overall_score", "r1_only"])
    if optimizer.log_output:
        print("Cleared branch history for {} site path".format(len(ins_file.get_all_sites())))

    # Recursively try to identify more candidate sites for mixing
    for leaf in initial.get_leaves():
        _do_site_mixing(leaf, tried.union(set(top_priority)), pairs, optimizer)


def try_anisotropy(initial, optimizer):
    """
    Test if adding anisotropy reduces R1 value.  If it does, do so.
    """

    new_ins = initial.get_res_copy()

    #  Try with anisotropy
    new_ins.add_anisotropy()
    annotation = "Added anisotropy"
    if optimizer.log_output:
        print("Trying step: {}".format(annotation))
    iteration = optimizer.history.run_iter(new_ins, initial, annotation)

    if iteration is not None:
        optimizer.history.save(iteration)
        if iteration.r1 > initial.r1:
            initial.propagate()


def try_exti(initial, optimizer):
    """
    Test if adding extinguishing reduces R1 value.  If it does, do so.
    """
    new_ins = initial.get_res_copy()

    #  Try with extinguishing
    new_ins.add_exti()
    annotation = "Added extinction"
    if optimizer.log_output:
        print("Trying step: {}".format(annotation))
    iteration = optimizer.history.run_iter(new_ins, initial, "Added extinction")

    # If exti helped, add it to the history
    if iteration is not None:
        optimizer.history.save(iteration)
        if iteration.r1 > initial.r1:
            initial.propagate()


def use_suggested_weights(initial, optimizer):
    """
    Stop re-initializing weights each time--use previously suggested weights
    """
    new_ins = initial.get_res_copy()

    #  Try with suggested
    new_ins.remove_command("WGHT")
    new_ins.commands.append(("WGHT", new_ins.suggested_weight_vals))
    annotation = "Used suggested weights"
    if optimizer.log_output:
        print("Trying step: {}".format(annotation))
    iteration = optimizer.history.run_iter(new_ins, initial, "Used suggested weights")

    if iteration is not None:
        optimizer.history.save(iteration)
        if iteration.r1 > initial.r1:
            initial.propagate()


def ensure_identified_elements_are_present(initial, optimizer):
    """
    Prune any branches that have missing elements
    """

    new_ins = initial.get_res_copy()

    identified_elements = set(map(lambda el: el.get_name(), new_ins.elements))
    result_elements = set()
    for site in new_ins.get_all_sites():
        result_elements.add(site.el_string.upper())
    missing = identified_elements.difference(result_elements)
    if len(missing) > 0:
        if optimizer.log_output:
            print("Truncated branch due to missing elements: {}".format(" ".join(missing)))
        optimizer.history.clean_history(n_to_keep=0, branch=initial)
