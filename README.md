| **About** | [![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE) [![Python](https://img.shields.io/badge/python-2.7%20-blue.svg)
| :------ | :------- |

# Satumut
`satumut` is a python package, wrappped around the [pmx package](https://github.com/deGrootLab/pmx) that performs saturated mutations on proteins to study their effects on protein-ligand interactions via [PELE simulations](http://www.nostrumbiodiscovery.com/pele.html).  

Given a position of a residue within a protein system:
1. it mutates to all the other 19 aminoacids by creating 19 PDbs
2. Then, it will create the files necessary for the PELE simulations the .yaml and the .sh files for each of the protein systems.

## Software requirements
* [numpy](https://numpy.org/)
* [scipy](https://www.scipy.org/)
* [matplotlib](https://matplotlib.org/)
* [pmx](https://pypi.org/project/pmx/)
* [seaborn](https://seaborn.pydata.org/)
* [fpdf](https://pyfpdf.readthedocs.io/en/latest/#installation)
* [biopython](https://biopython.org/wiki/Documentation)

