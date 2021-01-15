from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name="saturated_mutagenesis-etiur", author="Ruite Xiang", author_email="ruite.xiang@bsc.es",
      description="Study the effects of mutations on Protein-Ligand interactions",
      url="https://github.com/etiur/saturated_mutagenesis", license="MIT", version="0.0.1",
      packages=find_packages(), python_requires=">=2.7", long_description=long_description,
      long_description_content_type="test/markdown",
      classifiers=["Programming Language :: Python :: 2.7",
                   "License :: OSI Approved :: MIT License",
                   "Operating System :: Unix",
                   "Intended Audience :: Science/Research",
                   "Natural Language :: English",
                   "Environment :: Console",
                   "Development Status :: 1 - Planning",
                   ])