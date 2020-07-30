# illstack
Download:

pip install --user -e .


Running illstack:

Change basepath in example/istk-params_tng.txt for TNG simulation and example/istk-params_ill.txt for Illustris simulation to the directory storing the snapshot data.

Set base parameters in istk-params*.txt:

-Basepath for simulation data

-radial limits and number of bins for profiles

-searching for halos by mass option: Option 1 inputs upper/lower mass bounds given as mass_low and mass_high, Option 2 inputs a central mass and searches +- the mass_option_percent, given as a fraction

-scaled radius option, True for scaled (r/r200c) and False for unscaled (r)

-searching for halos by mass kind, options "stellar" or "halo" for stellar mass or halo mass


For use on multiple snapshots and/or mass ranges, use generate.py:

In the first few lines set the desired simulations (Illustris and TNG available) using sim (ill,tng), desired central masses, and desired population cut (for plotting) (high/low SFR, red/blue color, high/low stellar mass, or no cut) with cut (sfr,color,mstar,no_cut). 

It will automatically create shell scripts that call profiles_mult.py to create the npz files,then showstacks_*.py will plot.

Depending on individual permissions settings, the two subprocess.call in generate.py might give an error like "permission denied." Currently the way around this is to run generate.py once so the shell scripts are created even if they aren't run, give permission to run the shell scripts with

chmod +x *.sh

then run generate.py again
