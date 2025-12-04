# CLM-ml_v2.CHATS

This file describes some details for the multilayer canopy model
(CLM-ml.v2) as used in the CHATS simulations for April and May 2007.

This is the code used in:

Bonan, G. B., S. P. Burns, and E. G. Patton (2026). Beyond surface fluxes: Observational and computational
needs of multilayer canopy models – A walnut orchard test case. Agricultural and Forest Meteorology, 378,
110960, https://doi.org/10.1016/j.agrformet.2025.110960.

Supplied here are the source code and input datasets to run the simulations reported in the manuscript.
The 30-minute observation dataset created and used for model evaluation (see the Agric. For. Meteorol. manuscript
for documentation) is available separately at https://doi.org/10.5281/zenodo.17426258.

***Directory structure***

   1. The main multilayer canopy model code is in the directory:

      multilayer_canopy

      Model variables are defined in MLCanopyFluxesType.F90<br/>
      Model physical constants are defined in MLclm_varcon.F90<br/>
      Model parameters (array dimensions) are defined in MLclm_varpar.F90<br/>
      Model run control variables are defined in MLclm_varctl.F90<br/>
      The main driver code is MLCanopyFluxesMod.F90

      The code in this directory needs to be included with the CLM source code when coupled to CLM.

   3. There are two directories for simulations uncoupled to CLM:

      offline_driver     = driver code for the offline case. The main driver is CLMml.F90<br/>
      offline_executable = directory containing the makefile, executable, and namelists

      The offline code (uncoupled to CLM) has a namelist file read in
      offline_driver/controlMod.F90. This namelist file sets the tower site,
      the year and the month of the simulation, input files, and other options.

   4. The following directories are dummy CLM stub code only used in the offline case:

      clm_share           = dummy CLM stub code from: share/src<br/>
      clm_src_biogeophys  = dummy CLM stub code from: src/biogeophys<br/>
      clm_src_cpl         = dummy CLM stub code from: src/cpl/nuopc<br/>
      clm_src_main        = dummy CLM stub code from: src/main<br/>
      clm_src_utils       = dummy CLM stub code from: src/utils

      When coupled to CLM, the CLM code is used instead. Some CLM files need to be modified
      to couple with the multilayer canopy.

   5. The RSL psihat look-up file is in the directory:

      rsl_lookup_tables

***Run the model***

   Input files are provided to run the model for the CHATS tower site for April and May 2007. 
   A Makefile is provided to compile the code into the executable prgm.exe. The Makefile will have to be modified
   for the particular compiler. From offline_executable, use the command:

   ./prgm.exe < nl.CHATS7.05.2007

   or

   ./prgm.exe < nl.CHATS7.04.2007

   Input files are read from the input_files directory.
   Model output is written to the output_files directory.
   Example model output is provided (see subroutine output in offline_driver/CLMml_driver.F90 for file format).
