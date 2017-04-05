# Python library for automated single crystal refinement using Shelx

Authors: E. Antono and J. Ling at Citrine Informatics
Please contact: jling (at) citrine (dot) io with any questions.

This package is intended to be used with the SHELX code suite:
Sheldrick, George M. "A short history of SHELX." Acta Crystallographica Section A: Foundations of Crystallography 64.1 (2008): 112-122.

The purpose of this package is to automate many of the steps in single crystal refinement.  It
automatically detects split occupancy, deficiency, wrong element assignments, and anisotropy.

# Installation:

Clone this git repo:

git clone git@github.com:CitrineInformatics/crystal-refinement.git

Then, use:

python setup.py install

to install the package.

# Dependencies:

- Python packages: pymatgen, numpy

- SHELXTL xl.exe and xs.exe files (or xl, xs files for mac)
If these files are named differently in your distribution, please rename them to match.

- Citrination client:  https://github.com/CitrineInformatics/python-citrination-client
You will also need to create an account at citrination.com

# Example

TODO: Add an example/tutorial of use

