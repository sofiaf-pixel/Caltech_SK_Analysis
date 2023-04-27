# Caltech_SK_Analysis
This repository contains a set of scripts that are commonly used during detectors modules testing at Caltech.


To make load curves, PT plots, and Psat as a function of T with bolometer fit:

>> Run the bash script: "./find_G_all"

This script loops "sk_G_analysis.py" over rows and cols and biases.
If you want to fix beta value at 2, edit "find_G_all" with your favorite text editor and refer to "sk_G_analysis_beta2.py" (line 10).
However, be careful to save the previous plots in a separated folder before doing this, otherwise the script will overwrite them.

Note: you have to move the data in the same folder as the script.
- For that, just copy the full folder DATE, where DATE is the date at which you took data (yyyymmdd format)

Note: in "sk_G_analysis.py" (and "sk_G_analysis_beta2.py") script, you have to change:
- "RUN" to correspond to the current module number
- "DATE" to the date at which data has been taken (yyyymmdd format)

When you run ./find_G_all, a new "output" folder will be created in the main folder. Plots are saved in this folder.

--------------------------------------

To make the histograms:

>> Do "python make_histo.py"

This script is taking the occurrence for Gc, Tc, beta, RnTi, and Psat308 parameters to plot the histograms.

--------------------------------------

To make the tile maps:

>> Do "python plot_tilemap150.py"

This script in plotting the values of Gc, Tc, beta, RnTi, and Psat308 in a color scale for each detector of the module, according the Short Keck wiring diagram ("150ghz_mapping.xlsx" Excel file).

Note: these 2 scripts are looking for the .pkl files generated in the output folder. You must run "./find_G_all" first.

Note: in these 2 scripts, you have to change:
- "RUN" to correspond to the current module number
- "DATE" to the date at which data has been taken (yyyymmdd format)

Note: to plot histograms and tile maps for a fit with a beta fixed at 2, you have to run them again with the new .pkl files (after having run ./find_G_all refering to sk_G_analysis_beta2.py).
