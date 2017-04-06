"""
This script gives one example of how to call the optimizer.  It runs the optimization on the sample input
file input.ins and input.hkl.
"""

from optparse import OptionParser
from crystal_refinement.optimizer import Optimizer

# Parse command line arguments
parser = OptionParser()
parser.add_option("-l", "--path-to-xl",
                  help="Path to xl executable, including executable filename, e.g. path/to/file/xl.exe", dest="xl_path")
parser.add_option("-s", "--path-to-xs",
                  help="Path to xs executable, including executable filename, e.g. path/to/file/xs.exe", dest="xs_path")
parser.add_option("-w", "--use-wine", help="set to true if running Wine on Mac", default="False", dest="use_wine")
parser.add_option("-i", "--ins-path",
                  help= "Path to input ins file, without filename, e.g. /path/to/file/", dest="ins_path")
parser.add_option("-a", "--input-prefix", dest="input_prefix",
                  help="Prefix for input file. e.g. if input file is input.ins, then the prefix is 'input'")
parser.add_option("-z", "--output-prefix", dest="output_prefix",
                  help="Prefix for output file. This is where result of optimization will be written")

(options, args) = parser.parse_args()

# Define inputs to optimizer
xl_path = options.xl_path
xs_path = options.xs_path
use_wine = (options.use_wine == "True" or options.use_wine == "true")
ins_path = options.ins_path
input_prefix = options.input_prefix
output_prefix = options.output_prefix

# Create optimizer object and run on example
opt = Optimizer()
opt.run(path_to_xl=xl_path, path_to_xs=xs_path, ins_path=ins_path, input_prefix=input_prefix,
        output_prefix=output_prefix, use_wine=use_wine)


