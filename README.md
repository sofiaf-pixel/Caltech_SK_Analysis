# Caltech_SK_Analysis
This repository contains a set of scripts that are commonly used during detectors modules testing at Caltech.

--------------------------------------

First of all, move the data folder to the same folder as the scripts.
- For that, just copy the full folder DATE, where DATE is the date at which you took data (yyyymmdd format), to the folder where you save this set of scripts.

--------------------------------------

To do load curves, PT plots, and Psat as a function of T with bolometer fit:

>> Run the bash script: "./find_G_all"

If this doesn't work, run the command "chmod 777 find_G_all" to enable it as a bash script.
This script loops "sk_G_analysis.py" over rows and columns.

Note: In "find_G_all" script, you have to change:
- "RUN" to the number of the module that is tested
- "DATE" to the date at which data has been taken (YYYYMMDD format)

Note: Beta parameter is by default fixed to 2.
To fit the bolometer model with a variable value of beta, edit "sk_G_analysis.py" script with your favorite text editor and change "fit_beta" to "True" (line 28).

Note: The range of the load curves can be changed so that they all superpose. This corresponds to the region where detectors are superconducting.
For that, "fitrange["rnti_low"]" and "fitrange["rnti_hgh"]" must be modified (lines 115 and 116). These values are in non-physical units.

When you run ./find_G_all, a new "output" folder will be created in the main folder. Plots are saved in this folder.
If beta parameter isn't fixed to 2 (fit_beta=True), a folder named "Fit_beta" is generated to save the plots.

--------------------------------------

To plot the histograms:

>> Do "python make_histo.py"

This script is taking the occurrence for Gc, Tc, beta, RnTi, and Psat308 parameters to plot the histograms.

Note: This script is looking for the .pkl files generated in the output folder. You must run "./find_G_all" first.

Note: In "make_histo.py" script, you have to change:
- "RUN" to the name of the module that is tested
- "DATE" to the date at which data has been taken (YYYYMMDD format)

Note: To plot the histograms for a variable value of beta, just change "fit_beta" to "True" in "make_histo.py" script (line 23).
Before that, make sure you've run ./find_G_all for both options of fit_beta.

Note: The histogram's x-axis scale must represent all the values for Gc, Tc, beta, RnTi, and Psat308 but be as tight as possible.
For that, play with parameters in lines:
- 58 for Gc
- 81 for Tc
- 99 for beta
- 115 for RnTi
- 131 for Psat308
The first parameter in these lines is the beginning of the scale, the second is the end of the scale, and the third is the number of histograms that are plotted in this range (11 correspond to 10 histogram bars).
The parameters in the 2 following lines must follow the start and end of the scale set above.

The histogram plots are saved in output/DATE/RUN_Gplot/Histograms folder.
If fit_beta=True, Histograms folder is inside Fit_beta folder.

--------------------------------------

To plot the tile maps:

>> Do "python plot_tilemap150.py"

This script is plotting the values of Gc, Tc, beta, RnTi, and Psat308 in a color scale for each detector of the module, according to the Short Keck wiring diagram ("150ghz_mapping.xlsx" Excel file).

Note: This script is looking for the .pkl files generated in the output folder. You must run "./find_G_all" first.

Note: In "plot_tilemap150.py" script, you have to change:
- "RUN" to the name of the module that is tested
- "DATE" to the date at which data has been taken (YYYYMMDD format)

Note: To plot the tile maps for a variable value of beta, just change "fit_beta" to "True" in "plot_tilemap150.py" script (line 12).
Before that, make sure you've run ./find_G_all for both options of fit_beta.

Note: The tile maps color scale must be adjusted to the values of Gc, Tc, beta, RnTi, and Psat308.
For that, play with parameters in lines:
- 88/89 for Gc
- 77/78 for Tc
- 99/100 for beta
- 110/111 for RnTi
- 66/67 for Psat308
"vmin" parameter is the bottom of the color scale, and "vmax" is the top.

The tile maps are saved in output/DATE/RUN_Gplot/Tilemaps folder.
If fit_beta=True, tile maps folder is inside Fit_beta folder.
