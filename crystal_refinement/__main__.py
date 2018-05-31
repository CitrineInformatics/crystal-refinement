from __future__ import absolute_import

import yaml, ast
from argparse import ArgumentParser
from crystal_refinement.Optimizer import Optimizer

def parse_config(config_file):
    with open(config_file) as f:
        config = yaml.load(f)

        # SHELX Config section
        xl_path = config.get("path_to_xl")
        xs_path = config.get("path_to_xs")
        use_wine = config.get("use_wine", False)

        # Optimizer Config section
        r1_threshold = config.get("r1_threshold", 0.1)
        r1_similarity_threshold = config.get("r1_similarity_threshold", 0.0075)
        overall_score_ratio_threshold = config.get("overall_score_ratio_threshold", 1.3)
        occupancy_threshold = config.get("occupancy_threshold", 0.02)
        max_n_leaves = config.get("max_n_leaves", 50)

        # Bond Length Config section
        bond_lengths = config.get("bond_lengths", None)
        if bond_lengths is not None:
            bond_length_assert_message = "{} is an invalid configuration for bond_lengths. Please enter a list of tuples " \
                                         "in the form [('Ge', 'Mn', 2.5)]".format(bond_lengths)
            try:
                bond_lengths = ast.literal_eval(bond_lengths)
                assert type(bond_lengths) == list, bond_length_assert_message
                for bond in bond_lengths:
                    assert type(bond) == tuple \
                           and len(bond) == 3 \
                           and type(bond[0]) == str \
                           and type(bond[1]) == str \
                           and type(bond[2]) == float, bond_length_assert_message

            except ValueError:
                assert False, bond_length_assert_message
        use_ml_model = config.get("use_ml_model", False)
        citrination_api_key = config.get("citrination_api_key", None)

        # Mixing Pairs Config section
        mixing_pairs = config.get("mixing_pairs", None)
        if mixing_pairs is not None:
            mixing_pairs_assert_message = "{} is an invalid configuration for bond_lengths. Please enter a list of " \
                                          "tuples in the form [('Au', 'Os')]".format(bond_lengths)
            try:
                mixing_pairs = ast.literal_eval(mixing_pairs)
                assert type(mixing_pairs) == list, mixing_pairs_assert_message
                for mixing_pair in mixing_pairs:
                    assert type(mixing_pair) == tuple \
                           and len(mixing_pair) == 2 \
                           and type(mixing_pair[0]) == str \
                           and type(mixing_pair[1]) == str, mixing_pairs_assert_message

            except ValueError:
                assert False, mixing_pairs_assert_message

        # Score Config section
        score_weighting = config.get("score_weighting", 0.5)
        ensure_identified_elements = config.get("ensure_identified_elements", True)

        # Logging Config section
        n_results = config.get("n_results", 10)
        suppress_output = config.get("suppress_output", True)
        log_output = config.get("log_output", False)

        return {
            "path_to_xl": xl_path,
            "path_to_xs": xs_path,
            "use_wine": use_wine,
            "bond_lengths": bond_lengths,
            "mixing_pairs": mixing_pairs,
            "use_ml_model": use_ml_model,
            "citrination_api_key": citrination_api_key,
            "ensure_identified_elements": ensure_identified_elements,
            "r1_similarity_threshold": r1_similarity_threshold,
            "occupancy_threshold": occupancy_threshold,
            "r1_threshold": r1_threshold,
            "overall_score_ratio_threshold": overall_score_ratio_threshold,
            "score_weighting": score_weighting,
            "max_n_leaves": max_n_leaves,
            "n_results": n_results,
            "suppress_output": suppress_output,
            "log_output": log_output
        }


def main():
    """
    Command line interface for the optimizer
    """
    arg_parser = ArgumentParser()

    arg_parser.add_argument("-d", "--directory",
                            help="Path to the directory containing the initial ins file, without filename, "
                                 "e.g. /path/to/file/", dest="directory")
    arg_parser.add_argument("-i", "--input", dest="input_prefix",
                            help="Prefix for input file. e.g. if input file is input.ins, then the prefix is 'input'")
    arg_parser.add_argument("-o", "--output", dest="output_prefix",
                            help="Prefix for output file. This is where result of optimization will be written")
    arg_parser.add_argument("-c", "--config", dest="config_file",
                            help="YAML Configuration file specifying additional parameters for the optimizer")
    arg_parser.add_argument("-p", "--parse", dest="parse", action='store_true',
                            help="Print the values of the parameters read from the configuration file, including "
                                 "default values if parameters were not specified in the configuration file.")
    arg_parser.add_argument("-r", "--run", dest="run", action='store_true', help="Run the optimizer.")

    args = arg_parser.parse_args()
    if not (args.parse or args.run):
        print("Must specify either parse (-p) or run (-r).")
        quit()
    if args.config_file is None:
        print("No configuration file was specified")
        quit()

    config_file = args.config_file
    config = None
    if args.parse:
        config = parse_config(config_file)
        print("Optimizer configuration")
        for k, v in sorted(config.items(), key=lambda tup: tup[0]):
            print("{}={}".format(k, v))
    if args.run:
        input_directory = args.directory
        input_prefix = args.input_prefix
        output_prefix = args.output_prefix
        if input_directory is None or input_prefix is None or output_prefix is None:
            print("Input directory, input prefix and output prefix must all be specified when running the optimizer.")
            quit()
        if config is None:
            config = parse_config(config_file)
        # Create and run optimizer object
        opt = Optimizer(input_prefix=input_prefix,
                        output_prefix=output_prefix,
                        input_directory=input_directory,
                        **config)
        opt.run()

if __name__ == "__main__":
    main()