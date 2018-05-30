from setuptools import setup, find_packages

setup(
    name='crystal_refinement',
    version='0.1.0',
    packages=find_packages(),
    install_requires=['numpy',
                      'pymatgen',
                      'citrination-client>=4.0.0',
                      'graphviz'],
    url='citrine.io',
    license='BSD License',
    author='jling, eantono',
    author_email='erin (at) citrine (dot) io',
    description='Helps automate shelx for single crystal refinement',
    entry_points={
            'console_scripts': [
                  'crystal_refinement=crystal_refinement.__main__:main'
            ]
    }
)
