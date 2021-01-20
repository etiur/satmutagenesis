| **About** | [![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE) [![Python](https://img.shields.io/badge/python-2.7%20-blue.svg)
| :------ | :------- |

# Satumut
`satumut` is a python package, wrappped around the [pmx package](https://github.com/deGrootLab/pmx) that performs saturated mutations on proteins to study their effects on protein-ligand interactions via [PELE simulations](http://www.nostrumbiodiscovery.com/pele.html). It can also be used to design active sites (like a hydrolase site) by performing single mutagenesis around the bound substrate in multiple rounds.  

Given a position of a residue within a protein system:
1. It mutates to all the other 19 aminoacids by creating 19 PDBs
2. Then, it will create the files necessary for the PELE simulations the .yaml and the .sh files for each of the protein systems.

Given a PDB file of a protein-ligand system:
1. It mutates all the residues near the ester carbon in the substrate to Serine and selects the best after a PELE simulation.
2. A second round is made of the same approach to select the best positions for the Histidine of the catalytic triad.
3. A final round is made to select the best positions for the Aspartate/Glutamate of the catalytic triad.

## Software requirements
* [numpy](https://numpy.org/)
* [scipy](https://www.scipy.org/)
* [matplotlib](https://matplotlib.org/)
* [pmx](https://pypi.org/project/pmx/)
* [seaborn](https://seaborn.pydata.org/)
* [fpdf](https://pyfpdf.readthedocs.io/en/latest/#installation)

## Development

The scripts are continously modified and improved to have more functions and utilities.
