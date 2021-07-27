#!/usr/bin/env python

# imports
import os
import sys

# get the cwd were all the scripts are 
CWD = "/".join(__file__.split("/")[0:-1]); sys.path.insert(0, CWD)

# import the functions
import functions as fun

# general imports
import pandas as pd
import argparse, os
import pandas as pd
import numpy as np
from argparse import RawTextHelpFormatter
import shutil 

# paths under scripts
generate_plots  = "%s/generatePlots.R"%CWD

# get the paths of the colonyzer_fitness_env
CondaDir =  "/".join(sys.executable.split("/")[0:-4])
parametryzer = "%s/envs/colonyzer_fitness_env/bin/parametryzer"%CondaDir
colonyzer = "%s/envs/colonyzer_fitness_env/bin/colonyzer"%CondaDir

description = """
This is a pipeline to generate fitness analysis files and report from a set of images (each of them cotaining an image of a 96-well plate)contained in a specific directory (specified with the -i command). These images should have the first 15 characters as plate identifiers (error displayed if not) and the following should be a measurement of time, where if there is a "_" it indicates day_time (for example: ID_number_YYYYMMDD_HHMM, like img_0_20190821_2322, would be good).

This script should be run from an environment (colonyzer_fitness_env) created like this:

    conda create --name colonyzer_fitness_env python=2.7
    conda activate colonyzer_fitness_env
    pip install Colonyzer2
    pip install pygame
    conda install matplotlib
    # conda install -c anaconda pil --> not possible
    conda install -c anaconda pillow 
    conda install pandas
    conda install scipy
    pip install sobol
    #pip install opencv-python  --> took too long
    conda install -c conda-forge opencv
    #pip install --upgrade pillow

    Make sure that all pip runs are made from the proper environment

    However, it has to be run from colonyzer_fitness_env_run (also including the R packages), which needs:

    conda create --name colonyzer_fitness_env_run python=3.6
    conda activate colonyzer_fitness_env_run
    conda install pandas
    conda install -c anaconda biopython 
    conda install -c r r-sp 
    conda install -c conda-forge r-quantreg 
    conda install -c r r-jpeg
    conda install -c anaconda spyder 
    conda install matplotlib=3.3.0
    conda install seaborn
    conda install -c r r-knitr 
    conda install -c anaconda xlrd 

    R
    install.packages("DEoptim", repos="https://CRAN.R-project.org")
    install.packages("optparse")
    install.packages("qfa", repos="http://R-Forge.R-project.org")

    This is necessary because the program has to be run in python2.7 (colonyzer_fitness_env), but other things have to be run in python3 (colonyzer_fitness_env_run).

This program will try several combinations of image analysis provided by colonyzer. Each image (for 96-well plates) takes ~3 s to be analyzed.

The parameter combinations to be passed to colonyzer are:

- lc: Enable lighting correction
- diffims: If lighting correction switched on, attempt to correct for lighting differences between images in timecourse (can induce slight negative cell density estimates).
- edgemask: During lighting correction, use intensity gradient & morphology for final image segmentation instead of thresholding.
- cut: Cut culture signal from first image to make pseudo-empty plate
- greenlab: Check for presence of GreenLab lids on plates

It should be run on colonyzer_fitness_env_run.

"""

# arguments              
parser = argparse.ArgumentParser(description=description, formatter_class=RawTextHelpFormatter)

# mandatory args
parser.add_argument("-i", "--input", dest="input_dir", required=True, type=str, default=False, help="input dir, where all the images are")
parser.add_argument("-o", "--output", dest="output_dir", required=True, type=str, default=False, help="output dir, where all the files are written. Don't set it under inputdir, because the colonyzer pipeline recognizes as 'already analized' if there are any folders under input_dir that are output of the pipeline.")

# optional args
parser.add_argument("-n", "--n_wells", dest="n_wells", required=False, type=int, default=96, help="The number of wells of the plate")
parser.add_argument("-pcomb", "--parm_combinations", dest="parm_combinations", required=False, type=str, default="all", help="The parameter combinations to get passed to colonyzer. These can be any of lc,diffims,edgemask,cut,greenlab ; provided in a comma sepparated manner. If you provide 'all', it will run meaningful combinations of these.")
parser.add_argument("-man", "--manual_setting", dest="manual_setting", action="store_true", required=False, default=True, help="set manually the coordintes of the grid, which generates a Colonyzer.txt file with parametryzer")
parser.add_argument("--replace", dest="replace", action="store_true", required=False, default=False, help="Replace existing files")
parser.add_argument("-s", "--steps", dest="steps", required=False, type=str, default="image_correction,colonyzer,visualization,analysis", help="The steps to be run, comma sepparated. It can be any combination of image_correction,colonyzer,visualization")
parser.add_argument("--only_set_coordinates", dest="only_set_coordinates", action="store_true", required=False, default=False, help="Only set the coordinates of the plate (run parametryzer)")

opt = parser.parse_args()

# change the directory to the images file
os.chdir(opt.input_dir)

# make output if not done
fun.make_folder(opt.output_dir)

# get the steps to be analyzed
steps = set(opt.steps.split(","))


######### IMAGE CORRECTION ###########
if "image_correction" in steps: pass

######################################

######### COLONYZER RUN ###########
if "colonyzer" in steps:

    # initialize the commands to add to the colonyzer cmd
    extra_cmds = ""

    # run parametryzer to indicate the grid if manual option is indicated:
    if opt.manual_setting:

        colonizer_coordinates = "./Colonyzer.txt"
        if fun.file_is_empty(colonizer_coordinates) or opt.replace: fun.run_cmd_env(parametryzer, env="colonyzer_fitness_env")

        # add the command to initialize from coordinates
        extra_cmds+=" --initpos"

    # exit if only set coordinates
    if opt.only_set_coordinates: exit(0)

    # in the default manner, you have to infer meaningful combinations of parameters
    if opt.parm_combinations=="all":

        # define extra parameter combinations
        lc_parm_combinations = {("",), ("lc", "diffims", "edgemask"), ("lc",), ("lc", "diffims"), ("lc", "edgemask")}
        other_parm_combinations = {("",), ("cut", "greenlab"), ("cut",), ("greenlab",)}

        #all_parms = fun.make_flat_listOflists([set.union*([{tuple(set(lc_parms + other_parms).difference({""}))} for other_parms in other_parm_combinations]) for lc_parms in lc_parm_combinations])

        all_parms = [tuple(set(lc_parms + other_parms).difference({""})) for other_parms in other_parm_combinations for lc_parms in lc_parm_combinations]

    # if some are specified, get them
    else: all_parms = [tuple(opt.parm_combinations.split(","))]

    # define the image names that you expect
    image_names_withoutExtension = set([".".join(f.split(".")[0:-1]) for f in os.listdir(opt.input_dir) if f.split(".")[-1] in {"tif", "png", "jpg", "jpeg", "tiff"}])

    # go through these parms
    for parms in all_parms:

        # sort
        sorted_parms = sorted(parms)
        parms_str = "_".join(sorted_parms)
        extra_cmds_parmCombination = "".join([" --%s "%x for x in sorted_parms])

        # change the parms str if blank
        if parms_str=="": parms_str = "noExtraParms"

        # define the outdirs
        outdir = "%s/output_%s"%(opt.output_dir, parms_str)
        outdir_tmp = "%s_tmp"%outdir

        # check if all images have a data file (which means that they have been analyzed in outdir/Output_Data)
        all_images_analized = False
        Output_Data_dir = "%s/Output_Data"%outdir
        if os.path.isdir(Output_Data_dir):
            Output_Data_content = set(os.listdir(Output_Data_dir))

            if all(["%s.out"%f in Output_Data_content for f in image_names_withoutExtension]): all_images_analized = True

        # run the cmd with these parms
        if all_images_analized is False or opt.replace:
            print("\n\n\nwriting things for parameters %s"%parms_str)

            # create the folder (and empty if existing). Once you are here you will repeat the analysis of the data
            if os.path.isdir(outdir): shutil.rmtree(outdir)
            if os.path.isdir(outdir_tmp): shutil.rmtree(outdir_tmp)
            os.mkdir(outdir_tmp)

            # run colonizer, which will generate data under .
            colonyzer_cmd = "%s %s %s --plots --remove --fmt %i"%(colonyzer, extra_cmds, extra_cmds_parmCombination, opt.n_wells)
            fun.run_cmd_env(colonyzer_cmd, env="colonyzer_fitness_env")

            # once it is donde, move to outdir_tmp
            for folder in ["Output_Images", "Output_Data", "Output_Reports"]: fun.run_cmd("mv %s/ %s/"%(folder, outdir_tmp))

            # change the name, which marks that everything finished well
            os.rename(outdir_tmp, outdir)

###################################

########### VISUALIZATION #########
if "visualization" in steps: 

    # go through each of the directories of data and stack all .out files into one
    for outdir in os.listdir(opt.output_dir):
        if not os.path.isdir("%s/%s"%(opt.output_dir, outdir)): continue

        print("\n\nworking on %s"%outdir)

        # get the data path
        data_path = "%s/%s/Output_Data"%(opt.output_dir, outdir)

        # initialize a df
        all_df = pd.DataFrame()

        # add all the dfs of all images
        for f in [x for x in os.listdir(data_path) if x.endswith(".dat")]: 

            # get df
            df = pd.read_csv("%s/%s"%(data_path, f), sep="\t", header=None)

            # append
            all_df = all_df.append(df)

        # add barcode in the first place, instead of the filename
        all_df[0] = fun.get_barcode_for_filenames(all_df[0], filename_to_barcode_fn="initial_tests_Ewa")

        # sort the values
        all_df = all_df.sort_values(by=[0,1,2])

        # change the NaN by "NA"
        def change_NaN_to_str(cell):

            if pd.isna(cell): return "NA"
            else: return cell

        all_df = all_df.applymap(change_NaN_to_str)

        # write the csv under outdir
        all_df.to_csv("%s/%s/all_images_data.dat"%(opt.output_dir, outdir), sep="\t", index=False, header=False)

        # create the files that are necessary for the R qfa package to generate the output files

        # experiment descrption: file describing the inoculation times, library and plate number for unique plates. 
        exp_df = pd.DataFrame()

        # get all plates
        for I, plateBarcode in enumerate(set([x.split("-")[0] for x in all_df[0]])): 

            startTime = min(all_df[all_df[0].apply(lambda x: x.startswith(plateBarcode))][0].apply(lambda y: "-".join(y.split("-")[1:])))
            dict_data = {"Barcode":plateBarcode, "Start.Time":startTime, "Treatment":1, "Medium":"YPD" ,"Screen":"miki_screen", "Library":"CBS138", "Plate":I+1, "RepQuad":1}

            exp_df = exp_df.append(pd.DataFrame({k: {I+1 : v} for k, v in dict_data.items()}))

        # write
        exp_df.to_csv("%s/%s/ExptDescription.txt"%(opt.output_dir, outdir), sep="\t", index=False, header=True)

        # library description: where you state, for each plate (from 1, 2, 3 ... and as many plates defined in ExptDescription.Plate, the name and the ORF, if interestning)
        lib_df = pd.DataFrame()

        # define the rows and cols
        nWells_ro_NrowsNcols = {96:(8, 12)}

        for barcode, plateID in exp_df[["Barcode", "Plate"]].values:
            for row in range(1, nWells_ro_NrowsNcols[opt.n_wells][0]+1):
                for col in range(1, nWells_ro_NrowsNcols[opt.n_wells][1]+1):

                    # add to df
                    dict_data = {"Library":"CBS138", "ORF":"unk", "Plate":plateID, "Row":row, "Column":col, "Notes":""}
                    lib_df = lib_df.append(pd.DataFrame({k: {plateID : v} for k, v in dict_data.items()}))

        # write
        lib_df.to_csv("%s/%s/LibraryDescriptions.txt"%(opt.output_dir, outdir), sep="\t", index=False, header=True)

        # orf-to-gene
        pd.DataFrame(["unk", "unk_gene"]).transpose().to_csv("%s/%s/ORF2GENE.txt"%(opt.output_dir, outdir), sep="\t", index=False, header=False)

        # generate the plots with R
        print("Running R to generate plots")
        fun.run_cmd("%s %s"%(generate_plots, "%s/%s"%(opt.output_dir, outdir)))

if "analysis" in steps:

    # for several measures of fitness, plot a violin plot of the distribution of fitness, connecting each of the measurements

    # load the fit results into a df
    fit_df = pd.DataFrame()
    for outdir in os.listdir(opt.output_dir): 
        if not os.path.isdir("%s/%s"%(opt.output_dir, outdir)): continue

        # get df
        df = pd.read_csv("%s/%s/logRegression_fits.tbl"%(opt.output_dir, outdir), sep="\t")

        # add things
        df["type_calculation"] = ["_".join(outdir.split("_")[1:])]*len(df)
        df["sampleID"] = df.apply(lambda r: "%s_%s_%s"%(r["Barcode"], r["Row"], r["Column"]), axis=1)

        fit_df = fit_df.append(df)

    # write the fit_df for further usage
    fit_df.to_csv("%s/integrated_data.tbl"%opt.output_dir, sep="\t", index=False, header=True)

    # map the meaning of each GR estimate
    grEstimate_to_description_all = {"r": "Generalised logistic model rate parameter",
                                 "nAUC": "Numerical Area Under Curve. This is a model-free fitness estimate.",
                                 "nr": " Numerical estimate of intrinsic growth rate. Growth rate estimated by fitting smoothing function to log of data, calculating numerical slope estimate across range of data and selecting the maximum estimate (should occur during exponential phase)",
                                 "nr_t": "Time at which maximum slope of log observations occurs",
                                 "maxslp": "Numerical estimate of maximum slope of growth curve.",
                                 "maxslp_t": "Time at which maximum slope of observations occurs",
                                 "MDR": "Maximum Doubling Rate",
                                 "MDP": "Maximum Doubling Potential",
                                 "DT": "Doubling Time. Estimated from the fit parms at t0. May be biased if there is lag phase",
                                 "AUC": "Area Under Curve (from model fit)",
                                 "MDRMDP": "Addinall et al. style fitness",
                                 "DT_h": "max DT in hours. This is a numerical estimate from data",
                                 "rsquare": "rsquare between the fit and the data"}

    # define interesting measures to plot
    #interesting_GR_measures = {"r", "nAUC", "DT_h", "rsquare", "maxslp_t"}
    interesting_GR_measures = set(grEstimate_to_description_all)
    grEstimate_to_description  = {grEstimate : grEstimate_to_description_all[grEstimate] for grEstimate in interesting_GR_measures}

    for grEstimate, desc in grEstimate_to_description.items(): print(grEstimate, ": ",desc)

    # make plots
    fun.plot_growthRate_distributions(fit_df, grEstimate_to_description, "%s/comparativeAnalysisPlots"%opt.output_dir, plots={"correlation_mats", "violin_plots", "values_per_well"})

    # now everything that has lc on it
    #df_only_lc = fit_df[fit_df.type_calculation.apply(lambda x: "lc" in x)]
    #fun.plot_growthRate_distributions(df_only_lc, grEstimate_to_description, "%s/comparativeAnalysisPlots_onlyLC"%opt.output_dir, plots={"values_per_well", "violin_plots"})

    # everything with LC and no cut
    #df_only_lc_NOcut = fit_df[fit_df.type_calculation.apply(lambda x: "lc" in x and "cut" not in x)]
    #fun.plot_growthRate_distributions(df_only_lc_NOcut, grEstimate_to_description, "%s/comparativeAnalysisPlots_onlyLCandNOcut"%opt.output_dir, plots={"correlation_mats", "violin_plots"})
    # make subplots



###################################

"""
- lc: Enable lighting correction
- diffims: If lighting correction switched on, attempt to correct for lighting differences between images in timecourse (can induce slight negative cell density estimates).
- edgemask: During lighting correction, use intensity gradient & morphology for final image segmentation instead of thresholding.
- cut: Cut culture signal from first image to make pseudo-empty plate
- greenlab: Check for presence of GreenLab lids on plates
"""