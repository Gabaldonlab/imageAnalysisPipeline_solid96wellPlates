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

description = """
This script takes the several fitness measurements and calculates susceptibility measurements for different plate_sets.
It should be run on colonyzer_fitness_env_run.

"""

# arguments              
parser = argparse.ArgumentParser(description=description, formatter_class=RawTextHelpFormatter)

# mandatory args
parser.add_argument("--fitness_measurement_files", dest="fitness_measurement_files", required=True, type=str, default=False, help="A comma-sepparated string with the fitness measurement files of different plate sets")
parser.add_argument("--output_dir", dest="output_dir", required=True, type=str, default=False, help="The output dir")
parser.add_argument("--metadata_file", dest="metadata_file", required=True, type=str, default=False, help="a .tab file with the metadata")
parser.add_argument("--pseudocount_log2_concentration", dest="pseudocount_log2_concentration", required=False, type=float, default=0.1, help="A float that is used to pseudocount the concentrations.")
parser.add_argument("--min_nAUC_to_beConsideredGrowing", dest="min_nAUC_to_beConsideredGrowing", required=False, type=float, default=0.5, help="A float that indicates the minimum nAUC to be considered growing. This may depend on the experiment. This is added in the 'is_growing' field")
parser.add_argument("--min_points_to_calculate_resistance_auc", dest="min_points_to_calculate_resistance_auc", required=False, type=int, default=4, help="An integer number indicating the minimum number of points required to calculate the rAUC")


opt = parser.parse_args()

print("calculating susceptibility measurements for different plate sets")


########## PROCESS INPUTS ##########

# get the full paths
fitness_measurement_files_list = [fun.get_fullpath(x) for x in opt.fitness_measurement_files.split(",")]
output_dir = fun.get_fullpath(opt.output_dir)
metadata_file = fun.get_fullpath(opt.metadata_file)

fun.make_folder(output_dir)

# get the fitness measurements data
fitness_df = pd.DataFrame()
for plate_setI, fitness_file in enumerate(fitness_measurement_files_list):

	df = pd.read_csv(fitness_file, sep="\t")
	df["plate_set"] = plate_setI+1

	fitness_df = fitness_df.append(df).reset_index(drop=True)

# add the metadata
metadata_df = pd.read_csv(metadata_file, sep="\t")
joining_fields = ["plate_set", "plate", "row", "column"]
unwanted_fitness_df_fields = set(metadata_df.keys()).difference(set(joining_fields))
fitness_df = fitness_df[[c for c in fitness_df.keys() if c not in unwanted_fitness_df_fields]]
fitness_df = fitness_df.merge(metadata_df, on=joining_fields, validate="one_to_one", how="left")
fitness_df = fitness_df.sort_values(by=joining_fields)

# add fields
fitness_df["concentration"] = fitness_df.concentration.apply(float)
fitness_df["replicateID"] = "r" + fitness_df.row.apply(str) + "c" + fitness_df.column.apply(str)
fitness_df["sampleID"] = fitness_df.strain + "_" + fitness_df.replicateID
fitness_df["sampleID_and_condition"] = fitness_df.sampleID + "_" + fitness_df.condition
fitness_df["log2_concentration"] = np.log2(fitness_df.concentration + opt.pseudocount_log2_concentration)
fitness_df["is_growing"]  = fitness_df.nAUC>=opt.min_nAUC_to_beConsideredGrowing # the nAUC to be considered growing
fitness_df["bad_spot"] = fitness_df.bad_spot.apply(bool)

# debugs
print("checking that the inputs are OK")
if any(pd.isna(fitness_df.concentration)): raise ValueError("There can't be NaNs in concentration")
if any(pd.isna(fitness_df.condition)): raise ValueError("There can't be NaNs in condition")

for c in set(fitness_df.condition):
	df_c = fitness_df[fitness_df.condition==c]

	set_strainTuples = {tuple(df_c[df_c.concentration==conc].strain) for conc in set(df_c.concentration)}
	if len(set_strainTuples)!=1: raise ValueError("ERROR: This script expects the strains in each spot to be the same in all analyzed plates of the same condition (or drug). This did not happen for condition=%s. Check the provided --metadata_file."%c)
	
for c in set(fitness_df.condition):

	expected_nsamples = len(set(fitness_df[(fitness_df.condition==c)].sampleID))
	if sum((fitness_df.condition==c) & (fitness_df.concentration==0.0))!=expected_nsamples: raise ValueError("There should be %i wells with concentration==0 for condition==%s. Note that this script expects the strains in each spot to be the same in all analyzed plates of the same condition (or drug)."%(expected_nsamples, c))

####################################

######### GET THE DF WITH THE DRUG SUSCEPTIBILITY MEASUREMENTS ###########

# init variables
df_susceptibility = pd.DataFrame()
fitness_estimates  = ["K", "r", "nr", "maxslp", "MDP", "MDR", "MDRMDP", "DT", "AUC", "DT_h", "nAUC", "DT_h_goodR2"]

# get the fitness df with relative values (for each condition, the fitness relative to the concentration==0), and save these measurements
fitness_df = fun.get_fitness_df_with_relativeFitnessEstimates(fitness_df, fitness_estimates)
fitness_df.to_csv("%s/stacked_fitness_measurements.tab"%output_dir, sep="\t", header=True, index=False)

# get the susceptibility df for each sampleID
susceptibility_df = fun.get_susceptibility_df(fitness_df, fitness_estimates, opt.pseudocount_log2_concentration, opt.min_points_to_calculate_resistance_auc)
susceptibility_df.to_csv("%s/susceptibility_measurements.tab"%output_dir, sep="\t", header=True, index=False)

print("Congratulations!! The susceptibility data was correctly generated")
##########################################################################
