#!/usr/bin/env python

# imports
import os
import sys

# get the cwd were all the scripts are 
CWD = "/".join(__file__.split("/")[0:-1]); sys.path.insert(0, CWD)

# import the functions
import functions as fun

# define links
image_analysis_pipeline = fun.get_fullpath("%s/image_analysis_fitnesss_pipeline.py"%CWD)

description = """
This is a script that launches the image analysis pipeline for images were you have 4 plates for image. It creates under --input and --output 4 folders with the following names:

- plate1 (upper left)
- plate2 (upper right)
- plate3 (lower left)
- plate4 (lower right)

Each of the folders will contain info of the plate in each corner. Each of the folders in the output will be filled with soft links to the images and the locations for the corresponding quadrant

This script will also generate some analysis plots under --output if a metadata file is provided. It should be run on colonyzer_fitness_env_run.

"""

# general imports
import pandas as pd
import argparse, os
import pandas as pd
import numpy as np
from argparse import RawTextHelpFormatter
import shutil 

# arguments              
parser = argparse.ArgumentParser(description=description, formatter_class=RawTextHelpFormatter)

# mandatory args
parser.add_argument("-i", "--input", dest="input_dir", required=True, type=str, default=False, help="input dir, where all the images are. This script will create one folder for each pipeline here")
parser.add_argument("-o", "--output", dest="output_dir", required=True, type=str, default=False, help="output dir, where all the files are written. Don't set it under inputdir, because the colonyzer pipeline recognizes as 'already analized' if there are any folders under input_dir that are output of the pipeline. ")

# optional args
parser.add_argument("-n", "--n_wells", dest="n_wells", required=False, type=int, default=96, help="The number of wells of the plate")
parser.add_argument("-pcomb", "--parm_combinations", dest="parm_combinations", required=False, type=str, default="all", help="The parameter combinations to get passed to colonyzer. These can be any of lc,diffims,edgemask,cut,greenlab ; provided in a comma sepparated manner. If you provide 'all', it will run meaningful combinations of these.")
parser.add_argument("-man", "--manual_setting", dest="manual_setting", action="store_true", required=False, default=True, help="set manually the coordintes of the grid, which generates a Colonyzer.txt file with parametryzer")
parser.add_argument("--replace", dest="replace", action="store_true", required=False, default=False, help="Replace existing files")
parser.add_argument("-s", "--steps", dest="steps", required=False, type=str, default="image_correction,colonyzer,visualization,analysis,metadata_analysis", help="The steps to be run, comma sepparated. It can be any combination of image_correction,colonyzer,visualization")
parser.add_argument("-m", "--metadata_file", dest="metadata_file", required=False, type=str, default=None, help="A path to the metadata .tab table, which has plate, row, col ids as index and other columns as metadata. It should be a tab-sepparated file")


opt = parser.parse_args()

# set the extra args to the cmd
optional_args = "--n_wells %i --parm_combinations %s --steps %s"%(opt.n_wells, opt.parm_combinations, opt.steps)
if opt.manual_setting: optional_args += " --manual_setting"
if opt.replace: optional_args += " --replace"

# redefine the paths to be full paths
opt.input_dir = fun.get_fullpath(opt.input_dir)
opt.output_dir = fun.get_fullpath(opt.output_dir)
opt.metadata_file = fun.get_fullpath(opt.metadata_file)

# remove anything previously done
if opt.replace and os.path.isdir(opt.output_dir): shutil.rmtree(opt.output_dir)

# make output if not done
fun.make_folder(opt.output_dir)

# mappings
plateID_to_quadrantName = {1:"upper_left", 2:"upper_right", 3:"lower_left", 4:"lower_right"}

# define the image names
images_names = set([f for f in os.listdir(opt.input_dir) if f.split(".")[-1] in {"tif", "png", "jpg", "jpeg", "tiff"}])

# define a folder that will contain the softlinked_images
softlinked_input_dir = "%s/softlinked_images"%opt.output_dir; fun.make_folder(softlinked_input_dir)

# initialize the integrated table
integrated_table = "%s/integrated_data_allPlates.tbl"%opt.output_dir
if fun.file_is_empty(integrated_table):
#if True:
    # initialize a df that will contain all calculations
    all_df = pd.DataFrame()

    # iterate over the quadrants of the input dir
    for plateID in [1, 2, 3, 4]:

        print("Working on plate %i"%plateID)

        # define the dirs and make them
        plate_input_Dir = "%s/plate_%i"%(softlinked_input_dir, plateID) 
        plate_output_Dir = "%s/plate_%i"%(opt.output_dir, plateID) 
        for d in [plate_input_Dir, plate_output_Dir]: fun.make_folder(d)

        # soft link the images to the plate
        for image in images_names: 
            image_link = "%s/%s"%(plate_input_Dir, image)
            fun.soft_link_files("%s/%s"%(opt.input_dir, image), image_link)

        # generate the colonyzer cmd
        image_analysis_cmd = "%s -i %s -o %s %s"%(image_analysis_pipeline, plate_input_Dir, plate_output_Dir, optional_args)

        integrated_data_file = "%s/integrated_data.tbl"%(plate_output_Dir)

        if fun.file_is_empty(integrated_data_file):

            # run the colonyzer cmd, only to generate the colonyzer runs
            fun.run_cmd("%s --only_set_coordinates"%image_analysis_cmd)

            # run the full colonyzer cmd
            fun.run_cmd(image_analysis_cmd)


        if "analysis" not in opt.steps: continue

        # append the df
        df = pd.read_csv(integrated_data_file, sep="\t")
        df["plate"] =  [plateID]*len(df)
        all_df = all_df.append(df)


    if "analysis" not in opt.steps: 
        print("Warning. Exiting after the integration of data of the 4 plates. If you want to do this you should specify 'analysis' in --steps")
        sys.exit(0)

    # append the metadata
    if opt.metadata_file is not None: 
        metadata_df = pd.read_csv(opt.metadata_file, sep="\t")
        all_df = all_df.merge(metadata_df, how="left", left_on=["plate", "Row", "Column"], right_on=["plate", "row", "column"], validate="many_to_one")

    # correct the rsquare
    def get_rsquare_to0(rsq):

        if rsq>0: return rsq
        else: return 0.0

    all_df["rsquare"] = all_df.rsquare.apply(get_rsquare_to0)

    # get the correct DT_h
    maxDT_h = 25.0

    def get_DT_good_rsq(DT_h, rsq, rsq_tshd=0.95):

        if rsq>=rsq_tshd: return DT_h
        else: return maxDT_h

    all_df["DT_h_goodR2"] = all_df.apply(lambda r: get_DT_good_rsq(r["DT_h"], r["rsquare"]), axis=1)

    # write the integrated dataframe
    all_df.to_csv(integrated_table, sep="\t", index=False, header=True)

else: all_df = pd.read_csv(integrated_table, sep="\t")

if "metadata_analysis" not in opt.steps:
    print("skipping the analysis based on the provided metadata. You should add metadata_analysis into --steps if you want to do this.")
    sys.exit(0)

# if metadata is provided, make some plots
if opt.metadata_file is not None:

    # load df
    metadata_df = pd.read_csv(opt.metadata_file, sep="\t")

    # define and make plots dir
    PlotsDir = "%s/plots"%opt.output_dir
    if not os.path.isdir(PlotsDir): os.mkdir(PlotsDir)

    # define the metadata fields
    metadata_fields = [c for c in metadata_df.keys() if c not in {"plate", "row", "column"}]

    # define interesting measures to plot
    all_df["inv_DT_h_goodR2"] = 1 / all_df.DT_h_goodR2


    #interesting_GRestimates = ["r", "nAUC", "DT_h_goodR2", "rsquare", "AUC", "maxslp_t", "inv_DT_h_goodR2"]
    interesting_GRestimates = ["nAUC", "maxslp_t", "DT_h_goodR2", "rsquare", "inv_DT_h_goodR2"]

    # generate the plots
    #fun.generate_plots_many_plates(all_df, metadata_fields, fileprefix="%s/allData"%(PlotsDir), interesting_GRestimates=interesting_GRestimates)

    # only with lc
    #df_only_lc = all_df[all_df.type_calculation.apply(lambda x: "lc" in x)]
    #fun.generate_plots_many_plates(df_only_lc, metadata_fields, fileprefix="%s/onlyLC"%(PlotsDir), interesting_GRestimates=interesting_GRestimates)

    # everything with LC and no cut
    #df_only_lc_NOcut = all_df[all_df.type_calculation.apply(lambda x: "lc" in x and "cut" not in x)]
    #fun.generate_plots_many_plates(df_only_lc_NOcut, metadata_fields, fileprefix="%s/onlyLCnoCUT"%(PlotsDir), interesting_GRestimates=interesting_GRestimates)

    # only one of the parms:
    df_only_good = all_df[(all_df.type_calculation=="diffims_greenlab_lc") & (all_df.apply(lambda r: "blank" not in r, axis=1))]
    fun.generate_plots_many_plates(df_only_good, metadata_fields, fileprefix="%s/onlyDifLcGl"%(PlotsDir), interesting_GRestimates=interesting_GRestimates)

