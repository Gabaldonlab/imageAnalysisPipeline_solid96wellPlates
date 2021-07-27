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


            



    