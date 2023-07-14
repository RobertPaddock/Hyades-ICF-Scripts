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
This script creates the Lagrangian mesh over the capsule required for Hyades. The theory this is based on is described in the appendix of my thesis. Note that the meshing is a challenging problem, and the script may not perform well for capsules that are not close to those used in my thesis.

**File writer**



Batch analysis

**Hyades batch load**
Allows a batch run of files (produced using Hyades Batch Generator) to be analysed simultaneously. This file loads the data, and saves it in a file 'Results' in the folder where the data is stored, and also to a central repository location.

**Hyades batch plot**
Following hyades batch load, this allows the results from the previous batch run to be plotted. These can be plotted against the changing variables.

**AllData**
Reads in all the data stored in the central repository from hyades batch load, to form a large database of implosions. These can then be plotted against each other to identify e.g. highest gain implosions. Different criteria can be set e.g. only show sims with CR<16. Each sim is giving a unique identifier, so that the settings used in that sim can be accessed.



Other code:

**AuxHeatingResults_Git**
Shows how the Results.mat files generated from Hyades Batch Load for the auxiliary heating files can then be used to estimate gains, and produce comparisons of aux heating for different capsules



**Distinguishable colours**
Not created by me - downloaded from file exchange. Used in batch plots to generate plots with maximally distinguishable colours.

Copyright for distinguishable colours:

Copyright (c) 2010-2011, Tim Holy
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

