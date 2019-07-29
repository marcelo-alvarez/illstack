# illstack

Running illstack:

generate.py is the main overall script. 

In the first few lines of generate.py, set the desired simulations (Illustris and TNG available currently) using sim (ill,tng), desired central masses (with ranges of 20 percent above and below) with mass (1e12,1e13) and mcenter_power (12,13), and desired population cut (high/low SFR, red/blue color, high/low stellar mass) with cut (sfr,color,mstar,no_cut).

It will automatically create shell scripts that call profiles_mult.py to create the npz files,then showstacks_*.py will plot.

Depending on individual permissions settings, the two subprocess.call in generate.py might give an error like "permission denied." Currently the way around this is to run generate.py once so the shell scripts are created even if they aren't run, give permission to run the shell scripts with

chmod +x *.sh

then run generate.py again
