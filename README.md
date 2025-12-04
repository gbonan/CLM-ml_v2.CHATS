# CLM-ml_v2.CHATS

This file describes some details for the multilayer canopy model
(CLM-ml.v2) as used in the CHATS simulations for April and May 2007.

This is the code used in:

Bonan, G. B., S. P. Burns, and E. G. Patton (2026). Beyond surface fluxes: Observational and computational
needs of multilayer canopy models – A walnut orchard test case. Agricultural and Forest Meteorology, 378,
110960, https://doi.org/10.1016/j.agrformet.2025.110960.

Supplied are the source code and input datasets to run the simulations reported in the manuscript.
The 30-minute observation dataset created and used for model evaluation, as documented in the manuscript,
is available at https://doi.org/10.5281/zenodo.17426258.

1. Directory structure

   a. The main multilayer canopy model code is in the directory:

      1. multilayer_canopy

      i.   Model variables are defined in MLCanopyFluxesType.F90
      ii.  Model physical constants are defined in MLclm_varcon.F90
      iii. Model parameters (array dimensions) are defined in MLclm_varpar.F90
      iv.  Model run control variables are defined in MLclm_varctl.F90

      The main driver code is MLCanopyFluxesMod.F90

   b. There are two directories for simulations uncoupled to CLM:

      1. offline_driver     = driver code for the offline case. The main driver is CLMml.F90
      2. offline_executable = directory containing the makefile, executable, and namelists

      The offline code (uncoupled to CLM) has a namelist file read in
      offline_driver/controlMod.F90. This namelist file sets the tower site,
      the year and the month of the simulation, input files, and other options.

   c. The following directories are dummy CLM stub code only used in the offline case:

      1. clm_share           = dummy CLM stub code from: share/src
      2. clm_src_biogeophys  = dummy CLM stub code from: src/biogeophys
      3. clm_src_cpl         = dummy CLM stub code from: src/cpl/nuopc
      4. clm_src_main        = dummy CLM stub code from: src/main
      5. clm_src_utils       = dummy CLM stub code from: src/utils

      When coupled to CLM, the CLM code is used instead. Some CLM files need to be modified
      to couple with the multilayer canopy.

   d. The RSL psihat look-up file is in the directory:

      o rsl_lookup_tables

2. Run the model for a particular tower site.

   Input files are provided to run the model (offline) for CHATS tower site for April and May 2007.
   From offline_executable, use the command:

   ./prgm.exe < nl.CHATS7.05.2007

   or

   ./prgm.exe < nl.CHATS7.04.2007

   Input files are read from the input_files directory.
   Model output is written to the output_files_directory.
   Example model output is provided (see subroutine output in offline_driver/CLMml_driver.F90 for file format)
