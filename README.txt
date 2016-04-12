README

README file for ""
This .zip file contains R code and functions used for the landscape simulator in R and the source code for the implementation of FDIST2 in R.  The code was included "as is" at the time of publication.

LandscapeSimulatorCode Folder:

0Code_ALL_Select_C.R: R function for implementation of the demographies under neutral or selected evolution
Code_AllSourceFiles20130320: Source code for loading the C functions that were used for the recurrence equations and some other R functions that were used in the analysis.  This file also contains some commented code for installing packages that are necessary for analysis.
Code_MakeEnviCalls.R: Code for making the environments
Code_MakeRefugia.R: Code for making the starting conditions for the refugia.
Demographies: This folder contains the landscape matrices of starting populations sizes for the refugia models and the metadata files for the refugia.Environments:  This folder contains all the randomly generated environments. The file All_Environments_Metadata2.txt contains the information on which environments were used for each replicate set of simulations.MakeMetafile.R:  This file is R code for making a metadata file for a particular run.
MigMat.txt:  The migration/dispersal matrix used for the simulation.
SamplingSchemes/SchemeRandom1.txt:  The locations of the x-y coordinates used the sample the landscapes.
seed_X & seed_Y:  Seeds used for randomly choosing the sampling locations.
src:  This folder contains the C code for the implementation of the drift, migration, mutation, and selection recurrence equations.

--------------------------------------------------------------------------

srcFDIST2inR Folder:
changeExtension.R: An R function for changing the extension of a file.convertDatasets.R:  An R function for converting datasets.dependencies.R: Code for checking if packages need to be installed.
FODR_process.R: An example of calling the FDIST2 wrapper function.
FODRwrapper.R:  An example of code that goes through the individual steps of the FDIST2 analysis.FSTcalcs.R:  R functions to calculate various FST estimates that we compared at points in time.
FSTset.R: Function for setting the FST calculation to use for the analysis.
getAllelesDataset.R: Functions for rearranging allelic data.
getHeDataset.R: Functions for calculating expected heterozygosity.
getNumGens.R: Get number of generations to run based on IM parameters.
getPvalues.R: Functions used to get p-values and convert p-values into q-values.
getRunName.R: Function used to get a run name with a time-stamp.
getSampleSizeEachPop.R: Functions used to calculate sample sizes.IMsims.R: Functions used for simulating island model for neutral or selected loci.
makePlots.R: Functions used for making plots.myexport.R: Function for exporting.sampleVarianceFST.R: Function for calculating the sampling variance.sourceAll.R:  Load all source files.src/SimIM.c: C code for the island model implementationtictoc.R:  Functions for calculating the time spent on a process.
