# Hyades-ICF-meshing-and-analysis
Code for the meshing, production of input decks, and analysis of 3 layer ICF capsules

This is a version of the code I used to generate simulations for the large scale campaigns described in my thesis, and in the papers https://doi.org/10.1098/rsta.2020.0224 and https://doi.org/10.1017/S0022377822000265.

Some errors may exist, and the code will likely need to be adjusted for other purposes - but it should hopefully provide a useful start point.


In general:

**Hyades Batch Generator**
This script allows the user to specify the laser properties and capsule dimensions to be simulated. These must be of the laser profile form and capsule design described in the papers. The comments in the code should give an idea of how this can be done. Multiple files can be generated at any one time using the functionality in the script. This will call Analytic4Layers to do the meshing, and then produce the input decks in the MultiFile folder (this may need to be created and pointed to beforehand for the code to run).

**Hyades Front**
This script takes a single Hyades file, opens the data, and converts it into a file suitable for the following analysis file. The Analysis script is then called to analyse the data. The results are then shown (the Report script is also called to produce the large figure). Changing the 'figures' quantity determines how many plots are produced.

**Helios front**
As with Hyades front, but allows Helios simulations to be analysed instead.

**Hyades batch load**
Allows a batch run of files (produced using Hyades Batch Generator) to be analysed simultaneously. This file loads the data, and saves it in a file 'Results' in the folder where the data is stored, and also to a central repository location.

**Hyades batch plot**
Following hyades batch load, this allows the results from the previous batch run to be plotted. These can be plotted against the changing variables.

**AllData**
Reads in all the data stored in the central repository from hyades batch load, to form a large database of implosions. These can then be plotted against each other to identify e.g. highest gain implosions. Different criteria can be set e.g. only show sims with CR<16. Each sim is giving a unique identifier, so that the settings used in that sim can be accessed.



There are some variations of the Hyades Batch Generator file, which are designed for different purposes.

**Hyades batch generator two colour**
This file is designed to allow you to specify the wavelength of the main pulse sequence, and (if required) a second pulse sequence. This is the general file that should be used.

**Hyades batch generator aux heating**
This file is designed to allow you to easily explore the auxiliary heating of capsules.

