# imageAnalysisPipeline_solid96wellPlates

This is a pipeline to extract fitness measurements and plots for a set of images of 4 96-well solid plates with yeast colonies. It runs colonyzer and the R package qfa (https://research.ncl.ac.uk/qfa/#qfa). If you want to go deep on the data analysis it may be good for you to understand how qfa works from the provided link. The parameters used have been optimised for the GabaldonLab experimental setup.

## Installation

In order to use the pipeline you have to first install `conda` (or Anaconda) from https://docs.anaconda.com/anaconda/install/linux/#installation (follow instructions there). Anaconda is a software installer and manager of dependencies, which can also be very handy beyond this pipeline. We tested this on `conda` version 4.8.0, so that you may consider using the same version if something fails with the latest `conda` version.

Once you have `conda` installed, download this repository with:

`git clone https://github.com/Gabaldonlab/imageAnalysisPipeline_solid96wellPlates`

This will create a folder called `imageAnalysisPipeline_solid96wellPlates` with all the code of the pipeline and files necessary for installation. Move into this folder with:

`cd imageAnalysisPipeline_solid96wellPlates`

Beyonde the code in the repository (which is found in the folder `scripts`), this pipeline uses some R and python packages (dependencies). These can be installed automatically by creating two `conda` environments, called `colonyzer_fitness_env_run` and `colonyzer_fitness_env`. `conda` environments are like containers of software, which can keep dependencies within the environment isolated from other software in the computer. We built two environments for this pipeline because there are some dependencies that are incompatible with others.

You can create the environments with:

`conda env create --file ./installation/colonyzer_fitness_env.yml --name colonyzer_fitness_env`

`conda env create --file ./installation/colonyzer_fitness_env_run.yml --name colonyzer_fitness_env_run`

NOTE: These commands will remove any previous environment with the same name, so make sure that you don't have any such environments.

You have to install some extra dependencies which don't go by conda. First for the `colonyzer_fitness_env`:

`conda activate colonyzer_fitness_env`

`pip install Colonyzer2`

`pip install pygame`

`pip install sobol`

Then install some R packages in the `colonyzer_fitness_env_run`. First initialize R:

`conda activate colonyzer_fitness_env_run`

`R`

Then install some packages:

`install.packages("DEoptim", repos="https://CRAN.R-project.org")`

`install.packages("optparse")`

`install.packages("qfa", repos="http://R-Forge.R-project.org")`

Exit R (you don't need to save the workspace):

`q()`

If you got here it means that everything is ready to use.

## Running the pipeline on an example

This repository includes some data to test the pipeline as an example. These are a few images of Candida glabrata strains growing in four different media (everything is in the folder `testing`). You can run the pipeline on these images to check that the installation went well. The folder `testing_inputs` includes the inputs of the pipeline:

- A folder with the images (one image per timepoint), named as `_0_<YYYYMMDD>_<HHMM>.tif`. It is important that the images are named like this. For the example, this folder is in `testing/testing_inputs/images`.

- A tab-sepparated file with the metadata (information about what is in each well). This is a table where each row indicates a spot through a combination of the `plate` (a number between 1 and 4), `row` (a number between 1 and 8) and `column` (a number between 1 and 12). `strain` indicates what is in each spot, and `condition` shows the type of plate. For the example, this folder is in `testing/testing_inputs/metadata.tab`. You can create these files in excel and save as .csv, seeting the delimiter to tabs.

In order to run the pipeline you have to activate the `conda` environment with:

`conda activate colonyzer_fitness_env_run`

You can now run the analysis with:

`./scripts/run_image_analysis_4plates_per_image.py -i ./testing/testing_inputs/images --output ./testing/testing_output -pcomb greenlab,lc,diffims --steps colonyzer,visualization,analysis,metadata_analysis --metadata_file ./testing/testing_inputs/metadata.tab`

This pipeline is semiautomatic. It requires input from the user for defining the position of the upper-left and bottom-right wells of each plate. A window will appear showing the images for each of the four plates (plate 1: top-left, plate 2: top-right, plate 3: bottom-left, plate 4: bottom-right). You have to indicate the positions on all images and and save the coordinates as 96-well plate with 'g' when you finish. This is a bit of a pain and I hope they make this automatic in the future.

At the end this will create some files in `./testing/testing_output` with the results.

## Interpreting the output

There are two types of fitness estimates (found in the `integrated_data_allPlates.tbl` file) calculated from the data:

- Model-based fitness measurements: these are estimated by fitting a generalised logistic model (http://en.wikipedia.org/wiki/Generalised_logistic_function#Generalised_logistic_differential_equation) to the time-vs-cell density curve. The model parameters give us different fitness estimates. These estimates can be useful if we have some spots that did not reach stationary phase (to predict maximum cell density, for example) or we have mixed samples with different growth times. These don't work well if we have slow-growing spots or non-logistic curves (which may happen because there is cell death after reaching stationary phase). `K`, `r`, `g`, `v`, `MDR`, `MDP`, `DT`, `AUC`, `MDRMDP`, `rsquare` (see below) are related to such model fitting.


- Non parametric (or numeric) fitness measuremenets: these are calculated directly from the data, without assuming any underlying growth model. I generally use these (`nAUC` and `DT_h`) if we have experiments with the same growth times. `nAUC`, `nr`, `maxslp`, `maxslp_t`, `DT_h` and `DT_h_goodR2` (see below) are non-parametric measurements.

The pipeline outputs many files under directory specified with `--output`:

- `integrated_data_allPlates.tbl` is a tab-sepparated file with several fitness measurements for each well of the plate. Each row corresponds to a well, as indicated by the `plate`, `row` and `column` columns (note that these are also equivalent to the columns in the metadata file). These are the interesting fields:

	- `strain` and `condition` are those provided by user in the metadata file.
	
	- `Inoc.Time` is the time of innoculation formated as `YYYY-MM-DD_HH-MM-SS`. This is read from the image files.

	- `XOffset` and `YOffset` indicate the corrdinates of the well.

	- `K`, `r`, `g` and `v` are the parameters of a generalised logistic model that is fit to the data. You can check the qfa documentation if you want more precise information. `K` (maximum predicted cell density) and `r` (predicted growth rate) are fitness estimates that may be used.

	- `d0` is the normalised cell density of the first observation.

	- `nAUC` is the Numerical Area Under Curve. This is a model-free fitness estimate, directly calculated from the data.

	- `nr` is a numerical estimate of intrinsic growth rate. Growth rate estimated by fitting smoothing function to log of data, calculating numerical slope estimate across range of data and selecting the maximum estimate (should occur during exponential phase).

	- `maxslp` is a numerical estimate of maximum slope of growth curve, and `maxslp_t` is the time at which this maximum slope of observations occurs. `maxslp_t` is a way to calculate the lag phase.

	- `MDR` (Maximum Doubling Rate), `MDP` (Maximum Doubling Potential), `DT` (Doubling Time estimated from the model fit at t0, which may be biased if there is a lag phase), `AUC` (Area Under Curve) and `MDRMDP` (Addinall et al. style fitness) are several fitness estimates calculated from the model fit. You can check the qfa manual (http://qfa.r-forge.r-project.org/docs/qfa-manual.pdf) for more information

	- `rsquare` is the coefficient of determination (https://en.wikipedia.org/wiki/Coefficient_of_determination) between the model fit and the data. You can use it to determine which curves have a good model fit (i.e. rsquare > 0.95).

	- `DT_h` is a numerical estimate for the maximum doubling time, in hours. `DT_h_goodR2` is the same value but only for those spots with a good model fit (rsquare>0.95). For poorly fit curves the `DT_h_goodR2` is set to 25.0 (very high). This `DT_h_goodR2` can be used to have as non-growing the samples with weird curves.

- `plate_*/comparativeAnalysisPlots_0XXXXXXXXXX_valuesPerWell.pdf` is a heatmap plot with the fitness values per spot for different fitness estimates. There is one such plot per plate.

- `plate_*/comparativeAnalysisPlots_corrMat.pdf` is a heatmap with the correlation between the different fitness estimates across all wells. There is one such plot per plate.

- `plate_*/output_diffims_greenlab_lc/output_plots.pdf` has the growth curves with the model fitting. There is one such plot per plate.

- `plots/` includes several violinplots that describe the general trends in the 4 plates. The name of the plot indicates how is the violin plot organized. For example `plots/onlyDifLcGl_x=strain_hue=plate.pdf` is a plot where the x values are different strains and the colors (hue) indicate the plate number. These are automatic plots, which means that they won't be beautiful if there are a lot of different strains. They are just useful to scroll and have an overal description of the data. You can make any plots we have with the data in  `integrated_data_allPlates.tbl`.

## Get susceptibility measurements

This repository includes a script (`./scripts/measure_susceptibility_severalPlateSets.py`) to integrate the `integrated_data_allPlates.tbl` files of different experiments where the plates had a range of concentrations of a drug, and calculate susceptibility measurements for each strain. This script expects the strains in each spot to be the same in all analyzed plates of the same condition (or drug).

As an example, we'll run this script assuming that the testing data (in `./testing/testing_inputs/images`) is a set of plates with 4 anidulafungin (ANI) concentrations (0, 0.5, 1 and 2). In this example, we have only one plate set, but you could have any number of plate sets (for example, 2 plate sets with 8 concentrations in total). In order to run `measure_susceptibility_severalPlateSets.py` you need:

- A number of `integrated_data_allPlates.tbl` files (related to each of the plate sets that you want to analyse), one for each plate set. In our example we'll use the `./testing/testing_output/integrated_data_allPlates.tb`, generated in the section `Running the pipeline on an example`.

- A tab-sepparated file with the metadata (information about what is in each well). This is a table where each row indicates a spot through a combination of the `plate_set` (an number from 1 to the number of `integrated_data_allPlates.tbl` files provided), `plate` (a number between 1 and 4), `row` (a number between 1 and 8) and `column` (a number between 1 and 12). `strain` indicates what is in each spot, `condition` should be the drug of the plate and `concentration` would be the drug concentration. There has to be always a plate with `concentration`==0 for each `condition` to get the drug susceptibility measurements. For the example, this file is in `testing/testing_inputs/metadata_susceptibility.tab`. You'll see that the `concentration`==0 plate has also `condition`==ANI, which is how the script will know that the other ANI plates are related to this `concentration`==0. This will allow to run this script on many plate sets with different drugs assayed on different batches. For example, you could run this script on four plate sets (2 sets for ANI concentrations (first batch) and  2 sets for fluconazole (FLZ) concentrations (second batch)), so that each ANI or FLZ sample would have a different `concentration`==0 plate (the one corresponding to the same batch, as indicated by the `condition` field). Finally, the metadata df should also have a `bad_spot` field, which is a boolean (TRUE or FALSE) indicating whether that spot should be used for the susceptibility measurements (MIC and rAUC). As a general rule you can set it to FALSE, keeping all the wells. However, this can be very useful if you want the pipeline to remove some specific spots (i.e. due to contamination) of a plate keeping the measurements for the other spots. Sometimes, removing spots may make the susceptibility measurements impossible (i.e.: if you remove the highest concentration and there is no MIC in the lower concentrations you can't know the MIC). These 'unknown measurements' will appear as NaNs (which is fine because you can just ignore them in downstream analysis).

In order to run the script you have to activate the `conda` environment with:

`conda activate colonyzer_fitness_env_run`

You can now run the analysis with:

`./scripts/measure_susceptibility_severalPlateSets.py --fitness_measurement_files ./testing/testing_output/integrated_data_allPlates.tbl --output_dir ./testing/testing_output_susceptibility --metadata_file ./testing/testing_inputs/metadata_susceptibility.tab`

Note that in this case there is only one `integrated_data_allPlates` provided, which is related to `plate_set`==1 in `metadata_susceptibility.tab`. You can provide any number of plate sets. For example, for an experiment with 12 ANI concentrations (three plate sets), the command would look like:

`./scripts/measure_susceptibility_severalPlateSets.py --fitness_measurement_files <integrated_data_allPlates_set1>,<integrated_data_allPlates_set2>,<integrated_data_allPlates_set3> --output_file <output_file> --metadata_file <metadata_file>` 

Here the --metadata_file should have values 1,2,3 for the `plate_set` columns, and `ANI` in the `condition` column for all rows. The order of the `plate_set` should be the one provided with the comma-sepparated string in --fitness_measurement_files. 

This script will generate two files:

- `./testing/testing_output_susceptibility/stacked_fitness_measurements.tab` includes the stacked fitness measurement files (provided with --fitness_measurement_files), with some extra fields:

	- `replicateID` indicates the replicate as `r<row>c<column>`, which is useful to access one same sample across different plates. Derived from this there is the `sampleID` field, which includes `<strain>_<replicateID>`.

	- `log2_concentration` is the `log2(concentration + pseudocount)`. This pseudocount is default o 0.1, but it can be specified with the --pseudocount_log2_concentration argument (i.e. --pseudocount_log2_concentration 0.1).

	- `is_growing` is a boolean that indicates whether the spot is growing, which is necessary for the calculation of rAUC. Spots are considered to be growing if `nAUC` is above `min_nAUC_to_beConsideredGrowing` (default to 0.5, but it can be personalized with --min_nAUC_to_beConsideredGrowing 0.5, for example). Note that this --min_nAUC_to_beConsideredGrowing threshold may be defined by looking at the violinplots in `testing_output/plots/`.

	- `*_rel` columns include the relative fitness estimate measurement to the `concentration`==0 spot for several fitness estimates (`K`, `r`, `nr`, `maxslp`, `MDP`, `MDR`, `MDRMDP`, `DT`, `AUC`, `DT_h`, `nAUC` and `DT_h_goodR2`). These values can be useful to plot `concentration` vs fitness (i. e. `nAUC`) and see the inhibition curves.

- `./testing/testing_output_susceptibility/susceptibility_measurements.tab`. Each row includes several susceptibility measurements for each `condition` (or drug), `sampleID` (which is a combination of `strain` and `replicateID`) and `fitness_estimate`. These the important columns:

	- `fitness_estimate` indicates, for each row, which is the fitness estimate in which the MIC and rAUC values (see below) were calcualted. This table includes several RELATIVE estimates (`K`, `r`, `nr`, `maxslp`, `MDP`, `MDR`, `MDRMDP`, `DT`, `AUC`, `DT_h`, `nAUC` and `DT_h_goodR2`) and their corresponding ABSOLUTE estimates (`K_rel`, `r_rel`, `nr_rel`, `maxslp_rel`, `MDP_rel`, `MDR_rel`, `MDRMDP_rel`, `DT_rel`, `AUC_rel`, `DT_h_rel`, `nAUC_rel` and `DT_h_goodR2_rel`)

	- `MIC_*` is, for a given fitness_estimate, the minimum concentration at which you have a fitness<0.5, <0.75 or <0.90 (dependnig on the %). Note that this is only meaningful (in terms of drug susceptibility) for RELATIVE fitness estimates (those that are called like `*_rel`), where the fitness values are normalised to the `concentration`==0.

	- `rAUC_concentration` is the Area Under the concentration vs fitness Curve, divided by a 'maximum AUC' (where fitness==1 across all concentrations). Again, this is mostly meaningful for for RELATIVE fitness estimates (those that are called like `*_rel`), where the fitness values are normalised to the `concentration`==0. For those, it should be a value between 0 and 1.

	- `rAUC_log2_concentration` is similar to `rAUC_concentration`, but with the `log2_concentration` instead of `concentration`. As drug concentrations follow mostly a logarithmic range it can happen that the highest concentrations dominate the `rAUC_concentration` values. `rAUC_log2_concentration` gives more importance to all the concentartion range, and it is the one that we used in the C. glabrata paper of antifungal drug resistance.

	- `fitness_conc0` is the fitness at  `concentration`==0. It is only meaningful for ABSOLUTE fitness estimates (those that are NOT called like `*_rel`).


Below are some examples of which rows you can keep for different questions (for example, imagine that you have an R dataframe with the `susceptibility_measurements.tab` loaded called df_susceptibility):

- What is the MIC 50 for each sample in anidulafungin (ANI)? 

`df_susceptibility[df_susceptibility$fitness_estimate=="nAUC_rel" & df_susceptibility$condition=="ANI"]$MIC_50`

- What is the rAUC (resistance measurement) as we calculated in the C. glabrata resistance paper?

`df_susceptibility[df_susceptibility$fitness_estimate=="nAUC_rel" & df_susceptibility$condition=="ANI"]$rAUC_log2_concentration`

- What is the total biomass production (nAUC) in the absence of the drug?

`df_susceptibility[df_susceptibility$fitness_estimate=="nAUC" & df_susceptibility$condition=="ANI"]$fitness_conc0`


- What is the maximum doubling time in hours in the absence of the drug?

`df_susceptibility[df_susceptibility$fitness_estimate=="DT_h" & df_susceptibility$condition=="ANI"]$fitness_conc0`

These are some examples, but it is reccommended that you try different estimates and see which may fit better your experiment, and also beware that some filtering may be necessary (i.e. with the `bad_spot` field).


## Misc comments

These are some comments related to the development of this pipeline, only interesting for deverlopment:

- The conda environments where created with `conda env export --no-builds --from-history -n <env_name> --file <env_name>.yml` from a computer where these environments were generated and tested.