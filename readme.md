This is a readme file for supplemental data for "Computer-generated isotope model achieves experiment-level accuracy of fidelity for position-specific isotope analysis". This includes all the chemical mechanisms generated as well as the code used to convert those mechanisms into the figures in the paper. Here is a list of the folders each contains:

# input files

This folder contains the input file, `propane_rmg_input_file.py`, used by RMG to generate the full model. With RMG installed, the full isotope model can be run with the command `python $rmg/scripts/isotopes.py --kineticIsotopeEffect simple propane_rmg_input_file.py` and a mechanism without isotopes can be created with `python $rmg/rmg.py propane_rmg_input_file.py`. `$rmg` represents the path to the directory containing RMG-Py source code. 

# mechanisms

The propane and methane mechanisms used in this work are under the folder `mechanisms`. The file ending indicates the type of file. `inp` are CHEMKIN input files, `cti` are Cantera input files, `txt` are species dictionaries corresponding to the CHEMKIN input file for loading the CHEMKIN input file into RMG, and `csv` is a datatable comparing the cantera names to properties of isotopologue, and is provided so RMG installation is not necessary to perform analysis.

If the file name contains `no_isotopes`, then the file contains an RMG run without adding isotopes. If the file name contains `KIE`, then the file contains isotope labels with kinetic isotope effects. `no_KIE` indicates isotopes without kinetic isotope effects.

To use these mechanisms with the code provided, create four folders `three_reaction_model`, `six_reaction_model`, `drg_model`, and `full_model`, and place the corresponding isotopic mechanisms inside with the cantera file named `chem.cti`, the chemkin file called `chem.inp`, the species dictionary file as `species_dictionary.txt`, and the isotopologue datatable `isotopomer_cluster_info.csv`. 

# code

The folder `code` contains methods for converting between isotopologues and enrichment values and example scripts to help through the process.

`analysis_methods.py` contains methods for converting between isotopologues and enrichment values.

`cantera_tools.py` contains methods for running the cantera module. Comes from the github repository github.com/goldmanm/tools.

`prepare_model.py`  to create a file which maps isotopologues in cantera to specific enrichments, which is necessary for isotopic analysis. Running this file requires RMG to be installed. The output form this file is given in the mechanism sections. 

`create_paper_figures.py` is a script that simulates the models produced by RMG and creates all the figures in the main body of the paper. Code from this file provides an example of how simulations can be run and how to convert from enrichments to isotopologue fractions and vise versa. You must be in the `code` directory to run this script. Images are saved in the folder `results`.

## dependencies

The code in this section was tested with Python 2.7 and Python 3.6 with the following packages and their version numbers:

* numpy 1.11.3 (1.15.3 for Py3)
* pandas 1.0.2 (0.23.6 for Py3)
* cantera 2.3.0a3 (2.4.0 for Py3)
* seaborn 0.8 (0.9.0 for Py3)
* statsmodels 0.8.0 (0.9.0 for Py3)
* matplotlib 2.0.2 (2.2.3 for Py3)

All of these should packages except seaborn come with RMG installation. For the prepare_model.py function to work RMG must be installed. Documentation for this is available at reactionmechanismgenerator.github.io/RMG-Py/. All other analysis should be doable without RMG installation.

# exp_data

CSV files containing 'Gilbert' data from analyzing the figures and tables in Gilbert et al. Measurement of position-specific 13C isotopic composition of propane at the nanomole level. Geochemica et Cosmochemica Acta 177:205-216, 2016. The data originating from figures was approximated using [Engauge Digitizer](https://markummitchell.github.io/engauge-digitizer/). These are used in creating figures in the paper. 

# results

When you run `create_paper_figures.py`, it will create a results directory and save images here.


