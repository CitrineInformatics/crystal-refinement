from distutils.core import setup

setup(
    name='crystal_refinement',
    version='0.1.0',
    packages=['crystal_refinement'],
    install_requires=['numpy',
                      'pymatgen',
                      'citrination-client',
                      'graphviz'],
    url='citrine.io',
    license='BSD License',
    author='jling, eantono',
    author_email='erin (at) citrine (dot) io',
    description='Helps automate shelx for single crystal refinement'
)
