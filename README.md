# Hyades-ICF-Scripts

Here I have tried to share some of the code that I used to generate and analyse Hyades simulations for the large scale campaigns described in my thesis, and in the papers https://doi.org/10.1098/rsta.2020.0224 and https://doi.org/10.1017/S0022377822000265.
Some errors may exist, and the code will likely need to be adjusted for other purposes - but it should hopefully provide a useful start point.

The different scripts here call each other, so will all need to be saved in the same location. They have a number of features. I have tried to group similar scripts below:


Analysis:

**Hyades Front**
Allows the user to specify a single Hyades .cdf file. This is opened, and the data converted ready for analysis. The analysis script is then called.

**Helios Front**
As above, but for simulations performed in the altenative code Helios (the front scripts allow a shared analysis script to be used, despite the two codes).

**Analysis**
Takes the opened file, and runs analysis to plot the data, and calculate key variables (CR, gain, implosion velocity etc.). The adiabat reported by this script cannot be trusted. Calls the Report script for the final plots.

**Report**
Produces the large, full screen figure with a number of key plots.




File Generation:

**Hyades Batch Generator**
This script allows the user to specify the laser properties and capsule dimensions to be simulated. These must be of the laser profile form and capsule design described in the papers. The comments in the code should give an idea of how this can be done. Multiple files can be generated at any one time using the functionality in the script. This will call Analytic4Layers to do the meshing, and then produce the input decks in the MultiFile folder (this may need to be created and pointed to beforehand for the code to run).

**Meshing**

**File writer**



Batch analysis

**Hyades batch load**
Allows a batch run of files (produced using Hyades Batch Generator) to be analysed simultaneously. This file loads the data, and saves it in a file 'Results' in the folder where the data is stored, and also to a central repository location.

**Hyades batch plot**
Following hyades batch load, this allows the results from the previous batch run to be plotted. These can be plotted against the changing variables.

**AllData**
Reads in all the data stored in the central repository from hyades batch load, to form a large database of implosions. These can then be plotted against each other to identify e.g. highest gain implosions. Different criteria can be set e.g. only show sims with CR<16. Each sim is giving a unique identifier, so that the settings used in that sim can be accessed.

