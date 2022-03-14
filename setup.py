from setuptools import setup, find_packages
import satumut

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name="satumut", author="Ruite Xiang", author_email="ruite.xiang@bsc.es",
      description="Study the effects of mutations on Protein-Ligand interactions",
      url="https://github.com/etiur/satumut", license="MIT",
      version="%s" % satumut.__version__,
      packages=find_packages(), python_requires=">=3.7", long_description=long_description,
      long_description_content_type="text/markdown",
      classifiers=["Programming Language :: Python :: 3.7",
                   "Programming Language :: Python :: 3.8",
                   "License :: OSI Approved :: MIT License",
                   "Operating System :: Unix",
                   "Intended Audience :: Science/Research",
                   "Natural Language :: English",
                   "Environment :: Console",
                   "Development Status :: 1 - Planning",
                   "Topic :: Scientific/Engineering :: Bio-Informatics"],
      install_requires=["pmx @ git+https://github.com/cescgina/pmx.git@develop", "fpdf", "matplotlib", "numpy", "pandas", "seaborn", "biopython", "mdtraj"],
      keywords="protein engineering, bioinformatics, mutate proteins, simulations")
