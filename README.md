| **About** | [![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE) [![Python](https://img.shields.io/badge/python-2.7%20-blue.svg)] |
| :------ | :------- |

# Saturated_Mutagenesis
`Saturated_mutagenesis` is a python package, wrappped around the [pmx package](https://github.com/deGrootLab/pmx) that performs mutations on proteins to study their effects on protein-ligand interactions via [PELE simulations](http://www.nostrumbiodiscovery.com/pele.html).  
Given a position of a residue within a protein:
1. it mutates to all the other 19 aminoacids by creating 19 PDbs
2. it will create the files necessary for the PELE simulations, the .yaml and the .sh files.
