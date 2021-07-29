# These are the functions of the image analysis pipeline. It is intended to be run from the colonyzer_fitness_env_run environment in colonyzer_fitness_env_run

import os
import sys
import re
import string
import random
import pandas as pd
import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pickle
import itertools
from matplotlib import pyplot as plt
from sklearn.metrics import auc
import seaborn as sns
import itertools


def union_empty_sets(set_iterable):

    """It returns the union of an iterable of empty sets"""

    return set.union(*list(set_iterable) + [set()])

def run_cmd(cmd):

    """Runs a cmd and it stops if it did not work"""

    out_stat = os.system(cmd); 
    if out_stat!=0: raise ValueError("\n%s\n did not finish correctly. Out status: %i"%(cmd, out_stat))

def run_cmd_env(cmd, env="colonyzer_fitness_env"):

    """Runs cmd on a given environment"""

    # define the cmds
    CondaDir =  "/".join(sys.executable.split("/")[0:-4])
    SOURCE_CONDA_CMD = "source %s/etc/profile.d/conda.sh"%CondaDir
    cmd_prefix = "%s && conda activate %s &&"%(SOURCE_CONDA_CMD, env)

    # define the running
    run_cmd("bash -c '%s %s'"%(cmd_prefix, cmd))

def file_is_empty(path): 
    
    """ask if a file is empty or does not exist """
    
    if not os.path.isfile(path):
        return_val = True
    elif os.stat(path).st_size==0:
        return_val = True
    else:
        return_val = False
            
    return return_val

def make_flat_listOflists(LoL):

    return list(itertools.chain.from_iterable(LoL))


def make_folder(f):

    if not os.path.isdir(f): os.mkdir(f)


def get_fullpath(x):

    """Takes a path and substitutes it bu the full path"""

    # normal
    if x.startswith("/"): return x

    # a ./    
    elif x.startswith("./"): return "%s/%s"%(os.getcwd(), "/".join(x.split("/")[1:]))

    # others (including ../)
    else: return "%s/%s"%(os.getcwd(), x)

def remove_file(f):

    if os.path.isfile(f): 

        try: run_cmd("rm %s > /dev/null 2>&1"%f)
        except: pass


def soft_link_files(origin, target):

    """This function takes an origin file and makes it accessible through a link (target)"""

    if file_is_empty(target):

        # rename as full paths
        origin = get_fullpath(origin)
        target = get_fullpath(target)

        # check that the origin exists
        if file_is_empty(origin): raise ValueError("The origin %s should exist"%origin)

        # remove previous lisqnk
        try: run_cmd("rm %s > /dev/null 2>&1"%target)
        except: pass

        soft_linking_std = "%s.softlinking.std"%(target)
        print("softlinking. The std is in %s"%soft_linking_std)
        run_cmd("ln -s %s %s > %s 2>&1"%(origin, target, soft_linking_std))
        remove_file(soft_linking_std)

    # check that it worked
    if file_is_empty(target): raise ValueError("The target %s should exist"%target)

def filename_to_barcode_initial_tests_Ewa(filename):

    """example: img_0_2090716_1448"""

    # get plateID
    plateID = "".join(filename.split("_")[0:2])

    # get timestamp
    d = list(filename.split("_")[2])
    t = list(filename.split("_")[3])
    d_string = "%s%s%s%s-%s%s-%s%s"%(d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7])
    t_string = "%s%s-%s%s-00"%(t[0], t[1], t[2], t[3])

    return "%s-%s_%s"%(plateID, d_string, t_string)


def get_barcode_from_filename(filename, fn="initial_tests_Ewa"):

    """ It takes a filename and it generates a barcode that is plateID-YYYY-MM-DD_HH-MM-SS. The mapping is defined by fn"""

    fnName_to_fn  = {"initial_tests_Ewa": filename_to_barcode_initial_tests_Ewa}

    return fnName_to_fn[fn](filename)


def get_barcode_for_filenames(filenames_series, filename_to_barcode_fn="initial_tests_Ewa"):

    """Takes a series of filenames and passes them to get_barcode_from_filename to get barcoded values. The barcode cannot exceed 11 chars, so that it is stripped accoringly by this function"""

    # get barcoded items
    barcoded_names  = filenames_series.apply(lambda x: get_barcode_from_filename(x, fn=filename_to_barcode_fn))

    # get barcode
    unique_barcodes = set(barcoded_names.apply(lambda x: x.split("-")[0]))
    barcode_to_differentLengths = {b : len(b) for b in unique_barcodes if len(b)!=11}

    # initialize a map between the long and the short barcode
    oldBar_to_newBar = {x:x for x in unique_barcodes.difference(set(barcode_to_differentLengths))}

    # adapt the barcodes
    for barcode, len_barcode in barcode_to_differentLengths.items():

        # if larger
        if len_barcode>11: 

            # understand if there is a common suffix or a common prefix
            
            # common prefic
            if len(set([x[0:2] for x in unique_barcodes]))==1: newbarcode = barcode[len_barcode-11:]

            # common suffix
            elif len(set([x[-2:] for x in unique_barcodes]))==1: newbarcode = barcode[0:11]

            else: raise ValueError("No common suffix or prefix could be found. Please specify IDs that have 11 characters or less in the image filenames")

        # if smaller
        elif len_barcode<11: newbarcode = barcode + "X"*(11-len_barcode)

        # save
        oldBar_to_newBar[barcode] = newbarcode

    # check that everything was ok
    if len(oldBar_to_newBar)!=len(set(oldBar_to_newBar.values())): 
        print("The barcode transformation is:\n", oldBar_to_newBar)
        raise ValueError("The barcodes were not transformed correctly")

    # return
    return barcoded_names.apply(lambda x: "%s-%s"%(oldBar_to_newBar[x.split("-")[0]], "-".join(x.split("-")[1:])))


def plot_growthRate_distributions(df, GRestimate_to_description, fileprefix, group_var="type_calculation", plots={"violin_plots", "correlation_mats", "values_per_well"}):

    """Takes a dataframe (df) where each row is a sample and it has several growth rate calculations (as described in GRestimate_to_description). The df should have a grouping variable (group_var) which defines the different ways to calculate gr. It writes plots under fileprefix. sampleID indicates the samples that are the same across group_var."""


    # a plot where you have subplots. Each subplot for one field in GRestimate_to_description, it contains the distribution of values, and also the mappings between them
    if "violin_plots" in plots: 

        fig = plt.figure(figsize=(10, len(GRestimate_to_description)*3.5)); i=1

        # go through each GR estimate
        for GR, description in GRestimate_to_description.items():
            print("Plots for %s"%GR)

            # initialize subplot
            ax = plt.subplot(len(GRestimate_to_description), 1,i); i+=1;

            # add the lines that map each of the positions
            for sample in set(df.sampleID): plt.plot(group_var, GR, data=df[df.sampleID==sample], color="gray", linestyle="--", linewidth=0.2)

            # get the violin plots
            axvn = sns.violinplot(x=group_var, y=GR, data=df)

            # set the title
            ax.set_title(description[0:100])

            # delete the labels, only put in the last one
            ax.set_xticklabels([])

        # add the last xticklabels
        axvn.set_xticklabels(df[group_var].unique(), rotation=90)

        fig.tight_layout()  # otherwise the right y-label is slightly 
        fig.savefig("%s_violin.pdf"%fileprefix, bbox_inches='tight')
        plt.close(fig)

    # correlation matrix plot for each group_var, where you depict a correlation between all the different ways of looking at GR
    if "correlation_mats" in plots:

        # get all groups
        groups = set(df[group_var])

        fig = plt.figure(figsize=(20, 20)); i=1

        # go through them
        for group in groups:
            print("Working on %s"%group)

            # get a df only with the members of this group, and with the fields in GRestimate_to_description
            df_g = df[df[group_var]==group][list(GRestimate_to_description.keys())]

            # create a correlation mat 
            corr = df_g.corr(method="spearman")

            # Generate a mask for the upper triangle
            mask = np.zeros_like(corr, dtype=np.bool)
            mask[np.triu_indices_from(mask)] = True

            # Set up the subplot
            ax = plt.subplot(5, 4,i); i+=1;

            # Generate a custom diverging colormap
            cmap = sns.diverging_palette(220, 10, as_cmap=True)

            # Draw the heatmap with the mask and correct aspect ratio
            sns.heatmap(corr, mask=mask, cmap=cmap, vmax=1, center=0, square=True, linewidths=.5, cbar_kws={"shrink": .5})

            ax.set_title(group)
            ax.set_xticklabels(corr.keys(), rotation=90)

        fig.tight_layout()  # otherwise the right y-label is slightly 
        fig.savefig("%s_corrMat.pdf"%fileprefix, bbox_inches='tight')
        plt.close(fig)


    # value of each parameter at each well, one subplot for each of the values in GRestimate_to_description
    if "values_per_well" in plots:

        # go through each barcode, which should be a plate ID
        all_barcodes = set(df.Barcode)

        # get all groups
        groups = set(df[group_var])

        for barcode in all_barcodes:
            df_b = df[df.Barcode==barcode]

            fig = plt.figure(figsize=(len(GRestimate_to_description)*12*0.45, len(groups)*8*0.3)); i=1

            # go through all the groups 
            for i_g, group in enumerate(groups):
                print("Working on %s"%group)

                # get the df for this group
                df_g = df_b[df_b[group_var]==group]

                # go through them and plot subplots
                for GRestimate, description in GRestimate_to_description.items():
                    
                    # generate a df where you have, for GRestimate, the rowID as index and the colID as colname, and the values are the median value across all different "group_var"
                    col_to_row_to_GR = {}

                    
                    for col in set(df_b["Column"]):
                        for row in set(df_b["Row"]):
                            col_to_row_to_GR.setdefault(col, {}).setdefault(row, df_g[(df_g.Column==col) & (df_g.Row==row)][GRestimate].iloc[0])

                    df_gr = pd.DataFrame(col_to_row_to_GR)
                    
                    # generate plots
                    ax = plt.subplot(len(groups), len(GRestimate_to_description),i); i+=1;
                    
                    
                    cmap = sns.diverging_palette(220, 10, as_cmap=True)
                    sns.heatmap(df_gr, cmap=cmap, vmax=max(df_gr.max()), linewidths=.5, cbar_kws={"shrink": .5}, fmt=".2f", annot=True)
                    ax.set_title("%s-%s"%(group, GRestimate))    

            fig.tight_layout()  # otherwise the right y-label is slightly 
            fig.savefig("%s_%s_valuesPerWell.pdf"%(fileprefix, barcode), bbox_inches='tight')
            plt.close(fig)

def generate_plots_many_plates(df, metadata_fields, fileprefix, interesting_GRestimates, plots={"boxplot_grid"}):

    """Takes a df such as the integrated_data_allPlates.tbl and makes some plots. 
    
    - metadata_fields are columns in df which are metadata. In some plots, several combinations of these will be used 
    - interesting_GRestimates are the estimates that will be considered for the plot
    - plots indicates which plots to make

    """

    # make a df that contains only the interesting_GRestimates as a column, in a long format
    df = df[["plate", "row", "column", "type_calculation"] + metadata_fields + interesting_GRestimates]
    df_long = pd.melt(df, id_vars=(["plate", "row", "column", "type_calculation"] + metadata_fields), var_name='type_GRestimate', value_name='value')

    if "boxplot_grid" in plots:

        # a grid of jittered boxplots, where the rows are type_calculation and the cols are type_GRestimate. The X and hue can be any combination of plate and the metadata

        # define universal stuff
        all_types_calculations = list(set(df_long.type_calculation))
        all_GRestimates = list(set(df_long.type_GRestimate))

        # get combinations
        all_combinations = list(itertools.combinations(["plate"] + metadata_fields, 2))
        all_combinations += [(y,x) for x,y in all_combinations]

        # go through each combination of X and hue
        for X_hue_comb in all_combinations:

            print(X_hue_comb)

            x = X_hue_comb[0]; hue = X_hue_comb[1]
            filename = "%s_x=%s_hue=%s.pdf"%(fileprefix, x, hue)
            print("Working on %s"%filename)

            # calculate the number of boxes
            nboxes = len(set(df_long.apply(lambda r: (r[x], r[hue]), axis=1)))

            # initialize fig
            fig = plt.figure(figsize=(len(all_types_calculations)*nboxes*1, len(all_GRestimates)*3.5)); i=1
            #fig = plt.figure(); i=1

            # iterate through the columns

            for GRestimate in all_GRestimates:
                print(GRestimate)

                df_gr = df_long[df_long.type_GRestimate==GRestimate] 

                for typeCalc in all_types_calculations:
                    print(typeCalc)

                    df_calc = df_gr[df_gr.type_calculation==typeCalc]

                    # initialize subplot
                    ax = plt.subplot(len(all_GRestimates), len(all_types_calculations),i); i+=1

                    # add the boxplots plot
                    print("getting violin")
                    vn = sns.violinplot(x=x, y="value", hue=hue, data=df_calc)
                    #jit = sns.stripplot(x=x, y="value", hue=hue, data=df_calc, jitter=True, dodge=True, edgecolor='black', alpha=0.45, linewidth=1)
                    #jit = sns.swarmplot(x=x, y="value", hue=hue, data=df_calc, dodge=True, edgecolor='black', alpha=0.45, linewidth=1)

                    # avoid the jittering legend
                    #handles, labels = jit.get_legend_handles_labels()
                    #nhuevals = len(set(df_calc[hue]))
                    #jit.legend(handles[:nhuevals], labels[:nhuevals])

                    # subplot parms
                    ax.set_title("%s, %s"%(typeCalc, GRestimate))

            # make the violonplot
            #sns.catplot(x=x, y="value", hue=hue, data=df_long, kind="violin", row="type_calculation", col="type_GRestimate")
            #sns.violinplot(x=x, y="value", hue=hue, data=df_long)


            #plt.show()

            # save fig
            fig.tight_layout()  # otherwise the right y-label is slightly 
            fig.savefig(filename, bbox_inches='tight')
            plt.close(fig)

def get_fitness_df_with_relativeFitnessEstimates(fitness_df, fitness_estimates):

    """This function adds a set of *_rel fields to fitness_df, which are, for each condition, the fitness relative to the concentration==0 spot."""

    print("adding relative fitness to concentration==0")

    # correct the fitness estimates to avoid NaNs, 0s and infs
    fitEstimate_to_maxNoNInf = {fe : max(fitness_df[fitness_df[fe]!=np.inf][fe]) for fe in fitness_estimates}

    def get_correct_val(x, fitness_estimate):
        
        if x==np.inf: return fitEstimate_to_maxNoNInf[fitness_estimate]
        elif x!=np.nan and type(x)==float: return x
        else: raise ValueError("%s is not a valid %s"%(x, fitness_estimate))
        
    for fe in fitness_estimates: 
        
        # get the non nan vals
        fitness_df[fe] = fitness_df[fe].apply(lambda x: get_correct_val(x, fe))
        
        # add a pseudocount that is equivalent to the minimum, if there are any non negative values
        if any(fitness_df[fe]<0): 
            print("WARNING: There are some negative values in %s, modifying the data with a pseudocount"%fe)
            fitness_df[fe] = fitness_df[fe] + abs(min(fitness_df[fitness_df[fe]<0][fe]))    

    # define a df with the maximum growth rate (the one at concentration==0) for each combination of sampleID and assayed drugs
    df_max_gr = fitness_df[fitness_df.concentration==0.0].set_index("sampleID_and_condition", drop=False)[fitness_estimates]
    all_sampleID_and_condition = set(df_max_gr.index)

    fitEstimate_to_sampleIDandCondition_to_maxValue = {fe : {sampleIDandCondition : df_max_gr.loc[sampleIDandCondition, fe] for sampleIDandCondition in all_sampleID_and_condition} for fe in fitness_estimates}

    # add the relative fitness estimates
    fitness_estimates_rel = ["%s_rel"%x for x in fitness_estimates]

    def get_btw_0_and_1(x):
            
        if pd.isna(x): return 1.0
        elif x==np.inf: return 1.0
        elif x==-np.inf: return 0.0
        elif x<0: raise ValueError("there can't be any negative values")
        else: return x 

    fitness_df[fitness_estimates_rel] = fitness_df.apply(lambda r: pd.Series({"%s_rel"%fe : get_btw_0_and_1(np.divide(r[fe], fitEstimate_to_sampleIDandCondition_to_maxValue[fe][r["sampleID_and_condition"]])) for fe in fitness_estimates}), axis=1)

    return fitness_df


def get_MIC_for_EUCASTreplicate(df, fitness_estimate, concs_info, mic_fraction=0.5):

    """This function takes a df of one single eucast measurement, and returns the Minimal Inhibitory concentration, where the relative fitness is fitness_estimate, The df should be sorted by concentration."""

    # get the expected concs
    max_expected_conc = concs_info["max_conc"]
    first_concentration = concs_info["first_conc"]
    expected_conc_to_previous_conc = concs_info["conc_to_previous_conc"]

    # get the assayed concs
    assayed_concs = set(df.concentration)

    # calculate MIC
    concentrations_less_than_mic_fraction = df[df[fitness_estimate]<mic_fraction]["concentration"]

    # define a string for the warnings
    mic_string = "sampleID=%s|fitness_estimate=%s|MIC_%.2f"%(df.sampleID.iloc[0], fitness_estimate, mic_fraction)

    # define the real mic according to missing data

    # when there is no mic conc
    if len(concentrations_less_than_mic_fraction)==0:

        # when the max conc has been considered and no mic is found
        if max_expected_conc in assayed_concs: real_mic = (max_expected_conc*2)

        # else we can't know were the mic is
        else: 
            print("WARNING: There is no MIC, but the last concentration was not assayed for %s. MIC is set to NaN"%mic_string)
            real_mic = np.nan

    # when there is one
    else:

        # calculate the mic
        mic = concentrations_less_than_mic_fraction.min()

        # calculate the concentration before the mic
        df_conc_before_mic = df[df.concentration<mic]

        # when there is no such df, just keep the mic if there is only the first assayed concentration 
        if len(df_conc_before_mic)==0: 

            # if the mic is the first concentration or the first concentration is already below 0.5
            if mic==first_concentration: real_mic = mic    
            elif mic==0.0: real_mic = 0.01          
            else: 
                print("WARNING: We cound not find MIC for %s"%mic_string)
                real_mic = np.nan

        else:

            # get the known or expected concentrations
            conc_before_mic = df_conc_before_mic.iloc[-1].concentration
            expected_conc_before_mic = expected_conc_to_previous_conc[mic]

            # if the concentration before mic is not the expected one, just not consider
            if abs(conc_before_mic-expected_conc_before_mic)>=0.001: 
                print("WARNING: We cound not find MIC for %s"%mic_string)
                real_mic = np.nan
            else: real_mic = mic

    # if there is any missing 
    if real_mic==0: raise ValueError("mic can't be 0. Check how you calculate %s"%fitness_estimate)

    # debug
    """
    if pd.isna(real_mic):
        print(df)
        raise ValueError("MIC can't be nan")
    """

    return real_mic


def get_auc(x, y):

    """Takes an x and a y and returns the area under the curve"""

    return auc(x, y)


def get_AUC_for_EUCASTreplicate(df, fitness_estimate, concs_info, concentration_estimate, min_points_to_calculate_auc=4):

    """Takes a df were each row has info about a curve of a single EUCAST and returns the AUC with some corrections"""

    # get the expected concs
    max_expected_conc = concs_info["max_conc"]
    conc0 = concs_info["zero_conc"] 

    # calculate the auc if all the concentrations had a fitness of 1 (which is equal to no-drug if the fitness_estimate=1)
    max_auc = (max_expected_conc-conc0)*1

    # get the assayed concs
    assayed_concs = set(df[concentration_estimate])

    # when you lack less than 4 curves just drop
    if len(df)<min_points_to_calculate_auc: auc = np.nan

    # when they are all 0, just return 0
    elif sum(df[fitness_estimate]==0.0)==len(df): auc = 0.0

    else:

        # if the maximum growth is not measured, and the max measured is growing we discard the measurement, as it may change the results
        if max_expected_conc not in assayed_concs and df.iloc[-1].is_growing: auc = np.nan

        else:

            # get the values
            xvalues = df[concentration_estimate].values
            yvalues = df[fitness_estimate].values

            # if the initial is the first concentration, add the one
            if conc0 not in assayed_concs:
                xvalues = np.insert(xvalues, 0, conc0)
                yvalues = np.insert(yvalues, 0, 1.0)

            # calculate auc, relattive to the 
            auc = get_auc(xvalues, yvalues)/max_auc


    if auc<0.0: 

        print(df[[concentration_estimate, fitness_estimate, "is_growing"]])
        print(xvalues, yvalues)
        print(assayed_concs, conc0)
        raise ValueError("auc can't be 0. Check how you calculate %s"%fitness_estimate)

    return auc



def get_susceptibility_df(fitness_df, fitness_estimates, pseudocount_log2_concentration, min_points_to_calculate_auc):

    """Takes a fitness df and returns a df where each row is one sampleID-drug-fitness_estimate combination and there are susceptibility measurements (rAUC, MIC or initial fitness)"""

    # init the df that will contain the susceptibility estimates
    df_all = pd.DataFrame()

    # go through each condition
    for condition in sorted(set(fitness_df.condition)):
        print("getting susceptibility estimates for %s"%condition)

        # get the df
        fitness_df_c = fitness_df[fitness_df.condition==condition]

        # map each drug to the expected concentrations
        sorted_concentrations = sorted(set(fitness_df_c.concentration))
        concentrations_dict = {"max_conc":max(sorted_concentrations), "zero_conc":sorted_concentrations[0], "first_conc":sorted_concentrations[1], "conc_to_previous_conc":{c:sorted_concentrations[I-1] for I,c in enumerate(sorted_concentrations) if I>0}}

        sorted_log2_concentrations = [np.log2(c + pseudocount_log2_concentration) for c in sorted_concentrations]
        concentrations_dict_log2 = {"max_conc":max(sorted_log2_concentrations), "zero_conc":sorted_log2_concentrations[0], "first_conc":sorted_log2_concentrations[1], "conc_to_previous_conc":{c:sorted_log2_concentrations[I-1] for I,c in enumerate(sorted_log2_concentrations) if I>0}}

        # filter out bad spots
        fitness_df_c = fitness_df_c[~(fitness_df_c.bad_spot)]

        # go through all the fitness estimates (also the relative ones)
        for fitness_estimate in (fitness_estimates + ["%s_rel"%f for f in fitness_estimates]):
            print("fitness_estimate:", fitness_estimate)

            # define a grouped df, where each index is a unique sample ID
            grouped_df = fitness_df_c[["sampleID", "concentration", "is_growing", "log2_concentration", fitness_estimate]].sort_values    (by=["sampleID", "concentration"]).groupby("sampleID")

            # init a df with the MICs and AUCs for this concentration and fitness_estimate
            df_f = pd.DataFrame()

            # go through different MIC fractions
            for mic_fraction in [0.5, 0.75, 0.9]:

                df_f["MIC_%i"%(mic_fraction*100)] = grouped_df.apply(lambda x: get_MIC_for_EUCASTreplicate(x, fitness_estimate, concentrations_dict, mic_fraction=mic_fraction))

            # add the rAUC for log2 or not of the concentrations
            for conc_estimate, conc_info_dict in [("concentration", concentrations_dict), ("log2_concentration", concentrations_dict_log2)]:

                # get a series that map each sampleID to the AUC 
                df_f["rAUC_%s"%conc_estimate] = grouped_df.apply(lambda x: get_AUC_for_EUCASTreplicate(x, fitness_estimate, conc_info_dict, conc_estimate, min_points_to_calculate_auc=min_points_to_calculate_auc))

            # define the df for the initial fitness
            df_conc0 = fitness_df_c[(fitness_df_c["concentration"]==0.0)]
            df_f["fitness_conc0"] = df_conc0[["sampleID", fitness_estimate]].drop_duplicates().set_index("sampleID")[fitness_estimate]


            # keep df
            df_f = df_f.merge(fitness_df_c[["sampleID", "strain", "replicateID", "row", "column"]].set_index("sampleID").drop_duplicates(),  left_index=True, right_index=True, how="left",  validate="one_to_one")

            df_f["condition"] = condition
            df_f["fitness_estimate"] = fitness_estimate
            df_all = df_all.append(df_f)

    return df_all

