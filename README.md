# Python library for automated single crystal refinement with Shelx

Authors: E. Antono and J. Ling at Citrine Informatics

Please contact: jling (at) citrine (dot) io with any questions.

This package is intended to be used with the SHELX code suite:
Sheldrick, George M. "A short history of SHELX." Acta Crystallographica Section A: Foundations of Crystallography 64.1 (2008): 112-122.

The purpose of this package is to automate many of the steps in single crystal refinement.  It
automatically detects split occupancy, deficiency, wrong element assignments, and anisotropy.

# Installation:

Clone this git repo:
```sh
git clone git@github.com:CitrineInformatics/crystal-refinement.git
```
Then, use:
```sh
python setup.py install
```
to install the package.

# Dependencies:

- Python packages: pymatgen, numpy

- SHELXTL xl.exe and xs.exe files (or xl, xs files for mac).
If these files are named differently in your distribution, please rename them to match.

- Citrination client:  https://github.com/CitrineInformatics/python-citrination-client
You will also need to create an account at citrination.com

# Example

To run an example optimization on the ins file in the example folder:  
```sh
cd example  
python example.py --path-to-xl /path/to/xl.exe --path-to-xs /path/to/xs.exe --path-to-ins /path/to/insfile  
--input_prefix input --output_prefix output  
```

For example, I might run:  
```sh
python example.py --path-to-xl ./xl --path-to-xs ./xs --path-to-ins ../data/ --input-prefix input  
--output-prefix output  
```
This will generate an output file output.res containing the finalized crystal structure.

*Please note:* 

--path-to-xl and --path-to-xs should include the name of the executable itself (e.g. xl.exe).  Depending on the source of your SHELX executables, these might be called shelxl and shelxs, xl.exe and xs.exe, or xl and xs.  These should all be equivalent.

--path-to-ins should not include the file name of the ins file.  It should just include the path to the directory.

# How it works

The optimizer tries the following steps in order:
1) Identify crystal sites (try adding or removing Q peaks)
2) Switch elements at each crystal site
3) Try site mixing at each site
4) Try partial occupancy at each site
5) Try extinguishing
6) Try anisotropy 

At each step, the optimization takes into account both the fit R1 value as well as the bond lengths.
The idea bond lengths can be determined by 3 different methods:
1) User input via the bond_length argument to the Optimizer class
2) A machine learning model for bond length hosted at citrination.com.  To use this model, the user must have a free account on the public Citrination site.
3) Covalent radii of the elements

The optimizer will follow multiple paths to see which path yields the best final fit.  For example, if two different elements both yield similar R1 values at a given crystal site, the optimizer will try both options to see which give the best final fit.  This results in a branching tree structure.
The optimizer tree graph can be visualized using the generate_plot() method.


