import seaborn as sns
import pandas as pd
import statistics
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from typing import Union
from pathlib import Path
import gseapy as gp
import upsetplot as usp
from collections import Counter
import re

def read_out_top_target_genes(comparison_main_dir, form:Union["just_genes", "whole_dataframes"]="just_genes", 
modus: Union["intersection", "union", "occurrence_filter"]="union", affinitiy_filter=0, format:Union["list",  "set"]=list):
    """
    Early, in the whole work, created function for reading out the affinity values for TGs and their affinity values or only the TGs. Later not used for TG analysis as function had to be improved. 
    But this one still needed for the affinity distribution plot
    """
    
    dir_list = [path.path for path in os.scandir(comparison_main_dir) if path.is_dir()]
    target_gene_dict = {}
    for dir_path in dir_list:
        comparison_name, chromatin_status = get_involved_celltypes_and_chromatin_status(dir_path)
        if chromatin_status == "custom-HM-combination":
            continue
        if not comparison_name in target_gene_dict:
            target_gene_dict[comparison_name] = {}
        
        target_gene_dict = read_out_target_genes(comparison_path=dir_path, tf_tg_dict=target_gene_dict, comparison=comparison_name, form=form, mode=modus, affinitiy_filter=affinitiy_filter)
        
    if form == "just_genes":
        target_gene_df = pd.DataFrame.from_dict({(comparison, cell_type): target_gene_dict[comparison][cell_type] 
                                                                                for comparison in target_gene_dict.keys()
                                                                                for cell_type in target_gene_dict[comparison].keys()},
                                                                                orient="index")

    if form == "whole_dataframes":
        target_gene_df = pd.DataFrame.from_dict({(comparison, cell_type, hm): target_gene_dict[comparison][cell_type][hm] 
                                                                                for comparison in target_gene_dict.keys()
                                                                                for cell_type in target_gene_dict[comparison].keys()
                                                                                for hm in target_gene_dict[comparison][cell_type].keys()},
                                                                                orient="index")

    return target_gene_df

def read_out_target_genes(comparison_path: Union[str, Path], tf_tg_dict, comparison, form:Union["just_genes", "whole_dataframes"]="just_genes", mode:Union["intersection", "union", "occurrence_filter"]="union",  affinitiy_filter=0, format:Union["list",  "set"]=list):
    """Takes a path, checks it if it contains a directory with the top  TFs of the TF Prio run,
    and if so, reads the top TFs out and returns them as a list.
    """
    target_genes_dir = os.path.join(comparison_path, "output/output/org_exbio_tfprio_steps_distributionAnalysis_TopKTargetGenes/output/")
    if not os.path.isdir(target_genes_dir):
        print("The top target genes directory "+str(target_genes_dir)+" does not exist")
        return tf_tg_dict
    else:
        if not os.listdir(target_genes_dir):
            print("The top target genes directory "+str(target_genes_dir)+" is empty")
            return None
        else:
            for cell_type in os.listdir(target_genes_dir):
                cell_type_tg_dir = os.path.join(target_genes_dir, cell_type)
                if not cell_type in tf_tg_dict[comparison]:
                    tf_tg_dict[comparison][cell_type] = {}
                for hm in os.listdir(cell_type_tg_dir):
                    hm_cell_type_tg_dir = os.path.join(cell_type_tg_dir, hm)
                    if form == "just_genes":
                        tg_set = []
                    if form == "whole_dataframes":
                        tg_set= {}
                        tg_set["gene"] = []
                        tg_set["affinity"] = []
                    for tf_tg_file in os.listdir(hm_cell_type_tg_dir):
                        file_path = os.path.join(hm_cell_type_tg_dir, tf_tg_file)
                        if os.stat(file_path).st_size > 0:
                            file = pd.read_csv(file_path, sep="\t")
                            if form == "just_genes":
                                if affinitiy_filter > 0:
                                    file = file[file["Affinity"]>=affinitiy_filter]
                                new_tgs = file["Gene"].tolist()
                                if not tg_set:
                                    tg_set = new_tgs
                                if mode == "intersection":
                                    tg_set = tg_set.intersection(new_tgs)
                                
                                elif mode == "union" or mode == "occurrence_filter":
                                    tg_set.extend(new_tgs)

                            if form == "whole_dataframes":
                                tg_set["gene"].extend(file["Gene"])
                                tg_set["affinity"].extend(file["Affinity"])
                    if format == "set":
                        tg_set = set(tg_set)
                    if mode == "occurrence_filter":
                        gene_counts = col.Counter(tg_set)
                        tg_set = [gene for gene, count in gene_counts.items() if count > 1]

                    tf_tg_dict[comparison][cell_type][hm] = tg_set
    return tf_tg_dict


def get_involved_celltypes_and_chromatin_status(path):
    # improved version
    #e.g. path = "/nfs/data3/IHEC/TF_PRIO/comparisons/cells/b-cells-healthy-b-cells-leukemia-H3K27ac-H3k36me3-H3K4me1-H3K4me3-mRNA-totalRNA"
    # by splitting the folder name with the -"H3K27", its ensured that the part with the cell types of each comparison is used as a whole 
    # and no mixing up happens,  therefor  e.g. that a b-cells-healthy-monocyte-healthy would become monocyte-healthy-b-cells-healthy.
    # This is important for later functions as they iterate through the folders and then operate just on the folder 
    # which has the exact same cell types substring in it
    
    comparison = path.split("/")[-1] 
    comparison_cell_types = comparison.split("-H3K27")[0]
    
            
    if "H3K27ac-H3k36me3-H3K4me1-H3K4me3" in path or "H3K27ac-H3K36me3-H3K4me1-H3K4me3" in path:
        chromatin_status = "active_marks"
    elif "H3K27me3-H3K9me3" in path:
        chromatin_status = "suppressive_marks"
    else:
        #(maybe): TODO - improve custom HM combination detection
        print("The given path "+str(path)+" does no contain standard acitve or suppressive HM. \n Instead of the chromatin status, the combination of histon markers will be used as name.")
        chromatin_status = "custom-HM-combination"
        
    return comparison_cell_types, chromatin_status



def read_out_top_tfs_per_hm(comparison_main_dir):
    """As input directory a main comparisonfodler should be given,
    e.g. /nfs/data3/IHEC/TF_PRIO/comparisons/cells. Based on the input folder
    it iterates through every subdirectory and if this subdirectory contains the top TF
    list files (dcg.tsv), it will read them out and, together with the celltypes and chromatin status,
    and return all TFs as lists as value and with the according name of the folder as key in a dictionary.
    """
    
    dir_list = [path.path for path in os.scandir(comparison_main_dir) if path.is_dir()]
    overlap_df = pd.DataFrame(columns=['H3K4me1','H3K4me3','H3K27ac', 'H3K36me3', 'H3K27me3', 'H3K9me3'])
    for dir_path in dir_list:
        celltype1, celltype2, chromatin_status = get_involved_celltypes_and_chromatin_status(dir_path)
        if chromatin_status == "custom-HM-combination":
            continue
        comparison_name = "-".join([celltype1, celltype2])
        
        top_tf_comparison_name = "-".join([celltype1, celltype2, chromatin_status])
        top_tfs_dict = read_out_top_tf_hm_from_comparison_folder(dir_path)
        if top_tfs_dict is not None:
            if not comparison_name in overlap_df.index.tolist():
                overlap_df.loc[comparison_name] = [None, None, None, None, None, None]
            if chromatin_status == "suppressive_marks":
                top_tfs_dict.pop("H3K4me3", None)
            for hm_name,  top_tfs in top_tfs_dict.items():
                if top_tfs is not None:
                    overlap_df.at[comparison_name, hm_name] = set(top_tfs)
                else:
                    overlap_df.at[comparison_name, hm_name] = None

    return overlap_df

def read_out_top_tf_hm_from_comparison_folder(comparison_path: Union[str, Path]):
    #old version
    """Takes a path, checks it if it contains a directory with the top  TFs of the TF Prio run,
    and if so, reads the top TFs out and returns them as a list.
    """
    top_tf_hm_dict = {}
    #print(comparison_path)
    top_tf_hm_dir = os.path.join(comparison_path, "output/output/org_exbio_tfprio_steps_distributionAnalysis_CalculateDcgPerHm/output/")
    if not os.path.isdir(top_tf_hm_dir):
        print("The DcgRankPerHm directory "+str(top_tf_hm_dir)+" does not exist")
        return None
    else:
        if not os.listdir(top_tf_hm_dir):
            print("The DcgRankPerHm directory "+str(top_tf_hm_dir)+" is empty")
            return None
        else:
            for top_tf_hm_file in os.listdir(top_tf_hm_dir):
                hm_name = top_tf_hm_file.split(".")[0]
                top_tf_hm_file_path = os.path.join(top_tf_hm_dir, top_tf_hm_file)
                print(top_tf_hm_file_path)
                if os.stat(top_tf_hm_file_path).st_size > 0:
                    df = pd.read_csv(top_tf_hm_file_path, sep="\t", header=None)
                    top_tf_hm_dict[hm_name] = df[0]
                else:
                    top_tf_hm_dict[hm_name] = None
                
                
    return top_tf_hm_dict


def read_gtf_and_create_ensg_hgnc_dict(with_gene_id_version_number = False):
    path_to_gtf = "/nfs/data3/IHEC/TF_PRIO/additional_files/gencode.v29.annotation.gtf"
    with open(path_to_gtf) as f:
        gtf = list(f)
    gtf = [x for x in gtf if not x.startswith('#')]
    gtf_filtered = [x for x in gtf if 'gene_id "' in x and 'gene_name "' in x]
    if with_gene_id_version_number:
        gtf_ensg_hgnc_list = list(map(lambda x: (x.split('gene_id "')[1].split('"')[0], x.split('gene_name "')[1].split('"')[0]), gtf_filtered))
    else:
        gtf_ensg_hgnc_list = list(map(lambda x: (x.split('gene_id "')[1].split('.')[0], x.split('gene_name "')[1].split('"')[0]), gtf_filtered))
    gtf_ensg_hgnc_set = set(gtf_ensg_hgnc_list)
    gtf_dict = dict(gtf_ensg_hgnc_set)
    return gtf_dict


###### TODO: correct function
def create_overlap_df(comparison_main_dir):
    """As input directory a main comparisonfodler should be given,
    e.g. /nfs/data3/IHEC/TF_PRIO/comparisons/cells. Based on the input folder
    it iterates through every subdirectory and if this subdirectory contains the top TF
    list files (dcg.tsv), it will read them out and, together with the celltypes and chromatin status,
    and return all TFs as lists as value and with the according name of the folder as key in a dictionary.
    """
    
    dir_list = [path.path for path in os.scandir(comparison_main_dir) if path.is_dir()]
    overlap_df = pd.DataFrame(columns=['suppressive_marks','active_marks','overlap'])
    for dir_path in dir_list:
        celltype1, celltype2, chromatin_status = get_involved_celltypes_and_chromatin_status(dir_path)
        if chromatin_status == "custom-HM-combination":
            continue
        comparison_name = "-".join([celltype1, celltype2])
        
        top_tf_comparison_name = "-".join([celltype1, celltype2, chromatin_status])
        top_tfs = read_out_top_tf_from_comparison_folder(dir_path)
        if top_tfs is not None:
            if not comparison_name in overlap_df.index.tolist():
                overlap_df.loc[comparison_name] = [None, None, None]
            overlap_df.at[comparison_name, chromatin_status] = top_tfs

    return overlap_df


def read_out_top_tf_from_comparison_folder(comparison_path: Union[str, Path], modus: Union["list", "dataframe"]="list"):
    """Takes a path, checks it if it contains a directory with the top  TFs of the TF Prio run,
    and if so, reads the top TFs out and returns them as a list.
    """
    print(comparison_path)
    if not os.path.isdir(comparison_path):
        print("The used path "+str(comparison_path)+" does not exist")
        #Pass
        return None
    else:
        file_path = os.path.join(comparison_path,"output/output/org_exbio_tfprio_steps_distributionAnalysis_CalculateDcgRank/output/dcg.tsv")
        if not os.path.isfile(file_path):
            print("There is no dcg.tsv file for the given folder"+str(os.path.join(comparison_path,"output/output/org_exbio_tfprio_steps_distributionAnalysis_CalculateDcgRank/output")))
            #pass
            return None
        else:
            #print(file_path)
            if os.path.getsize(file_path) == 0:
                print("File at "+str(file_path)+" is empty. File will be deleted.")
                os.remove(file_path)
            else:
                print("Found top_tf file at: "+str(file_path))
                top_tf_df = pd.read_csv(file_path, sep="\t")
                if modus == "list":
                    return list(top_tf_df["TF"])
                elif modus == "dataframe":
                    return top_tf_df

def adjust_tfprio_used_TF(tf_list):
    """Takes the list of all used TFs by TFPrio and changes the names of the TFs in it.
    E.g. NFIC::TLX1 to NFIC and TLX1 or JUND(MA0492.1) to JUND. After this also removes duplicates.
    TODO: just as 
    """
    #new_tf_list = list(set([x.split("(")[0] for x in[x.split(":")[0] for x in tf_list]])) #old command
    #new_tf_list = list(set([y.split("(")[0] for x in tf_list for y in x.split("::")]))
    tf_without_parentheses = [re.sub(r'\(.*?\)', '', s) for s in tf_list]
    new_tf_list = list(set([y.split("(")[0] for x in tf_without_parentheses for y in x.split("::")]))
    return new_tf_list


def read_out_tpm_filtered_genes(comparison_dir, gtf_dict, tpm_filter_value=1, as_hgnc=False, as_united_list=False):
    tpm_valid_genes_dict = {}
    tpm_dir = os.path.join(comparison_dir,"output/output/org_exbio_tfprio_steps_rnaSeq_CalculateTPM/output/")
    if not os.path.isdir(tpm_dir):
        print("The tpm directory "+str(tpm_dir)+" does not exist")
        return tpm_valid_genes_dict
    else:
        if not os.listdir(tpm_dir):
            print("The tpm directory "+str(tpm_dir)+" is empty")
            return tpm_valid_genes_dict
        else:
            for cell_type_file in os.listdir(tpm_dir):
                cell_type = cell_type_file.split(".")[0]
                file_path = os.path.join(tpm_dir, cell_type_file)
                file = pd.read_csv(file_path, sep="\t")
                tpm_filtered_genes = file[file["mean"]>=tpm_filter_value]["gene_id"].tolist()
                if as_hgnc:
                    if as_united_list:
                        tpm_filtered_genes = list(convert_ensg_to_hgnc(tpm_filtered_genes, gtf_dict))
                    else:
                        tpm_filtered_genes = convert_ensg_to_hgnc(tpm_filtered_genes, gtf_dict)
                tpm_valid_genes_dict[cell_type]=tpm_filtered_genes
    if as_united_list:
        tpm_valid_genes_list = sum(tpm_valid_genes_dict.values(),[])
        return tpm_valid_genes_list
    return tpm_valid_genes_dict


def convert_ensg_to_hgnc(string_set, ensg_hgnc_dict):
    return {ensg_hgnc_dict.get(x) for x in string_set}


def read_out_filtered_target_genes(tf_list, comparison_path, modus:Union["all_genes", "top_k_genes"]="top_k_genes", affinity_filter_value=0.05):
    """Reads out the target genes above the set affinity_filter_value. Etther using the top_k_genes from the output folder or the files  from the input folder in 
    org_exbio_tfprio_steps_distributionAnalysis_TopKTargetGenes
    """
    list_of_filtered_genes = []
    base_dir = os.path.join(comparison_path, "output/output/org_exbio_tfprio_steps_distributionAnalysis_TopKTargetGenes")
    if modus == "top_k_genes":
        tf_files_list = [tf+".tsv" for tf in tf_list]
        #print(tf_files_list)
        gene_base_dir = os.path.join(base_dir, "output")
        for root, dirs, files in os.walk(gene_base_dir):
            for tf_file in files:
                if tf_file in tf_files_list:
                    #print(tf_file)
                    tf_file_path = os.path.join(root, tf_file)
                    file = pd.read_csv(tf_file_path, sep="\t")
                    file = file[file["Affinity"]>=affinity_filter_value]
                    new_tgs = file["Gene"].tolist()
                    list_of_filtered_genes.extend(new_tgs)

    elif modus == "all_genes":
        gene_base_dir = os.path.join(base_dir, "input")
        for root, dirs, files in os.walk(gene_base_dir):
            for tg_matrix in files:
                if tg_matrix == "dcg.tsv":
                    continue
                tg_matrix_path = os.path.join(root, tg_matrix)
                tg_matrix_df = pd.read_csv(tg_matrix_path,  sep="\t")
                for tf in tf_list:
                    if tf in tg_matrix_df.columns:
                        new_tgs = tg_matrix_df[tg_matrix_df[tf]>=affinity_filter_value]["Gene"].tolist()
                        list_of_filtered_genes.extend(new_tgs)

    return list_of_filtered_genes



def read_out_filtered_target_genes_dict(tf_df, all_comparisons_path, modus:Union["all_genes", "top_k_genes"]="top_k_genes", affinity_filter_value=0.05, tf_per_hm=False, per_cell_type=False):
    """Reads out the target genes above the set affinity_filter_value. Either using the top_k_genes from the output folder or the files  from the input folder in 
    org_exbio_tfprio_steps_distributionAnalysis_TopKTargetGenes

    """
    tf_dict = tf_df.to_dict("index")
    tf_filtered_tg_dict={}
    dir_list = [path.path for path in os.scandir(all_comparisons_path) if path.is_dir()]
    for comparison, values in tf_dict.items():
        if comparison not in tf_filtered_tg_dict.keys():
            tf_filtered_tg_dict[comparison] = {}
            for col_name in tf_df.columns:
                tf_filtered_tg_dict[comparison][col_name] = list()

        for dir in dir_list:
            if comparison in dir:
                comparison_path = os.path.join(all_comparisons_path, dir)
                base_dir = os.path.join(comparison_path, "output/output/org_exbio_tfprio_steps_distributionAnalysis_TopKTargetGenes")
                if modus == "top_k_genes":
                    for mode, tf_list in values.items():
                        tf_files_list = [tf+".tsv" for tf in tf_list]
                        if per_cell_type:
                            gene_base_dir = os.path.join(base_dir, "output", comparison)
                        else:
                            gene_base_dir = os.path.join(base_dir, "output")
                        for root, dirs, files in os.walk(gene_base_dir):
                            if mode == "active_marks" and any(hm_marker in dirs for hm_marker in ["H3K9me3", "H3K27me3"]):
                                    #print(root)
                                    #print("no TGs for: "+mode+" in dir: "+' '.join(dirs))
                                    continue
                            if mode == "suppressive_marks" and any(hm_marker in dirs for hm_marker in ["H3K4me1", "H3K4me3", "H3K27ac", "H3K36me3"]):
                                    #print(root)
                                    #print("no TGs for: "+mode+" in dir: "+' '.join(dirs))
                                    continue
                            for tf_file in files:
                                if tf_file in tf_files_list:
                                    tf_file_path = os.path.join(root, tf_file)
                                    file = pd.read_csv(tf_file_path, sep="\t")
                                    file = file[file["Affinity"]>=affinity_filter_value]
                                    new_tgs = file["Gene"].tolist()
                                    tf_filtered_tg_dict[comparison][mode].extend(new_tgs)
            
                elif modus == "all_genes":
                    if per_cell_type:
                        gene_base_dir = os.path.join(base_dir, "input", comparison)
                    else:
                        gene_base_dir = os.path.join(base_dir, "input")
                    for root, dirs, files in os.walk(gene_base_dir):
                        for tg_matrix in files:
                            if tg_matrix == "dcg.tsv":
                                continue
                            tg_matrix_path = os.path.join(root, tg_matrix)
                            tg_matrix_df = pd.read_csv(tg_matrix_path,  sep="\t")
                            for mode, tf_list in values.items():
                                if mode == "active_marks" and any(hm_marker in tg_matrix for hm_marker in ["H3K9me3", "H3K27me3"]):
                                    #print("no TGs for: "+mode+" in file: "+tg_matrix)
                                    continue
                                if mode == "suppressive_marks" and any(hm_marker in tg_matrix for hm_marker in ["H3K4me1", "H3K4me3", "H3K27ac", "H3K36me3"]):
                                    #print("no TGs for: "+mode+" in file: "+tg_matrix)
                                    continue
                                for tf in tf_list:
                                    if tf in tg_matrix_df.columns:
                                        new_tgs = tg_matrix_df[tg_matrix_df[tf]>=affinity_filter_value]["Gene"].tolist()
                                        tf_filtered_tg_dict[comparison][mode].extend(new_tgs)
                            
    return tf_filtered_tg_dict


def read_out_filtered_target_genes_df(tf_df, all_comparisons_path, modus:Union["all_genes", "top_k_genes"]="top_k_genes", affinity_filter_value=0.05, tf_per_hm=False, per_cell_type=False):
    """Reads out the target genes above the set affinity_filter_value. Either using the top_k_genes from the output folder or the files  from the input folder in 
    org_exbio_tfprio_steps_distributionAnalysis_TopKTargetGenes

    """
    tf_dict = tf_df.to_dict("index")
    tf_filtered_tg_dict={}
    dir_list = [path.path for path in os.scandir(all_comparisons_path) if path.is_dir()]
    for comparison, values in tf_dict.items():
        if comparison not in tf_filtered_tg_dict.keys():
            tf_filtered_tg_dict[comparison] = {}
            for col_name in tf_df.columns:
                tf_filtered_tg_dict[comparison][col_name] = list()

        for dir in dir_list:
            if comparison in dir:
                comparison_path = os.path.join(all_comparisons_path, dir)
                base_dir = os.path.join(comparison_path, "output/output/org_exbio_tfprio_steps_distributionAnalysis_TopKTargetGenes")
                if modus == "top_k_genes":
                    for mode, tf_list in values.items():
                        tf_files_list = [tf+".tsv" for tf in tf_list]
                        if per_cell_type:
                            gene_base_dir = os.path.join(base_dir, "output", comparison)
                        else:
                            gene_base_dir = os.path.join(base_dir, "output")
                        
                        for root, dirs, files in os.walk(gene_base_dir):
                            if tf_per_hm:
                                if os.path.basename(root) != mode:
                                    continue
                            else:
                                #if mode == "active_marks" and any(hm_marker in dirs for hm_marker in ["H3K9me3", "H3K27me3"]):
                                #        continue
                                #elif mode == "suppressive_marks" and any(hm_marker in dirs for hm_marker in ["H3K4me1", "H3K4me3", "H3K27ac", "H3K36me3"]):
                                #        continue          
                                if ((mode == "active_marks" and os.path.basename(root) in ["H3K4me1", "H3K4me3", "H3K27ac", "H3K36me3"]) | 
                                    (mode == "suppressive_marks" and os.path.basename(root) in ["H3K9me3", "H3K27me3"]) | (mode in ["unique", "intersection"])):    
                                    # prior change: for tf_file loop und alles danach ein tab zurück sodass auf  gleicher Ebene wie elif
                                    for tf_file in files:
                                        if tf_file in tf_files_list:
                                            tf_file_path = os.path.join(root, tf_file)
                                            file = pd.read_csv(tf_file_path, sep="\t")
                                            file = file[file["Affinity"]>=affinity_filter_value]
                                            new_tgs = file["Gene"].tolist()
                                            tf_filtered_tg_dict[comparison][mode].extend(new_tgs)
            
                elif modus == "all_genes":
                    if per_cell_type:
                        gene_base_dir = os.path.join(base_dir, "input", comparison)
                    else:
                        gene_base_dir = os.path.join(base_dir, "input")
                    for root, dirs, files in os.walk(gene_base_dir):
                        for tg_matrix in files:
                            if tg_matrix == "dcg.tsv":
                                continue
                            tg_matrix_path = os.path.join(root, tg_matrix)
                            tg_matrix_df = pd.read_csv(tg_matrix_path,  sep="\t")
                            for mode, tf_list in values.items():
                                if tf_per_hm:
                                    if mode not in tg_matrix:
                                        continue
                                else:
                                    if mode == "active_marks" and any(hm_marker in tg_matrix for hm_marker in ["H3K9me3", "H3K27me3"]):
                                        continue
                                    elif mode == "suppressive_marks" and any(hm_marker in tg_matrix for hm_marker in ["H3K4me1", "H3K4me3", "H3K27ac", "H3K36me3"]):
                                        continue
                                for tf in tf_list:
                                    if tf in tg_matrix_df.columns:
                                        new_tgs = tg_matrix_df[tg_matrix_df[tf]>=affinity_filter_value]["Gene"].tolist()
                                        tf_filtered_tg_dict[comparison][mode].extend(new_tgs)
    tg_df =  pd.DataFrame.from_dict(tf_filtered_tg_dict, orient="index")        
               
    return tg_df


def add_overlap_union_unique_shared_with_all_or_once(df):
    tf_overlap_df = df
    tf_overlap_df["intersection"] = tf_overlap_df.apply(lambda row: set(set(row["suppressive_marks"])&set(row["active_marks"])), axis=1)
    tf_overlap_df["union"] = tf_overlap_df.apply(lambda row: set(set(row["suppressive_marks"])|set(row["active_marks"])), axis=1)

    def find_unique_values(row):
        unique_values = []
        for value in row["union"]:
            if value not in row["intersection"]:
                if all(value not in other_list for other_list in tf_overlap_df["union"] if other_list != row["union"]):
                    unique_values.append(value)
        return set(unique_values)

    tf_overlap_df["unique"] = tf_overlap_df.apply(find_unique_values, axis=1)


    def find_shared_values(row):
        other_lists = set().union(*tf_overlap_df.loc[tf_overlap_df.index != row.name, 'union'])
        return set([value for value in row['union'] if value in other_lists])

    tf_overlap_df['shared_at_least_once'] = tf_overlap_df.apply(find_shared_values, axis=1)

    common_values = set.intersection(*map(set, tf_overlap_df['intersection']))

    def filter_common_values(row):
        return set([value for value in row['union'] if value in common_values])

    tf_overlap_df['shared_with_all'] = tf_overlap_df.apply(filter_common_values, axis=1)

    return tf_overlap_df


def read_out_top_tf_from_comparison_folder_with_tpm(comparison_path: Union[str, Path], gtf_dict,tpm_filter=0, modus: Union["list", "dataframe", "set"]="set"):
    """Takes a path, checks it if it contains a directory with the top  TFs of the TF Prio run,
    and if so, reads the top TFs out and returns them as a list.
    """
    tpm_filtered_valid_genes_hgnc = []
    tpm_filtered_valid_genes_hgnc = read_out_tpm_filtered_genes(comparison_dir=comparison_path, gtf_dict=gtf_dict, tpm_filter_value=tpm_filter, as_hgnc=True, as_united_list=True)
    if not os.path.isdir(comparison_path):
        print("The used path "+str(comparison_path)+" does not exist")
        return None
    else:
        file_path = os.path.join(comparison_path,"output/output/org_exbio_tfprio_steps_distributionAnalysis_CalculateDcgRank/output/dcg.tsv")
        if not os.path.isfile(file_path):
            print("There is no dcg.tsv file for the given folder"+str(os.path.join(comparison_path,"output/output/org_exbio_tfprio_steps_distributionAnalysis_CalculateDcgRank/output")))
            return None
        else:
            if os.path.getsize(file_path) == 0:
                print("File at "+str(file_path)+" is empty. File will be deleted.")
            else:
                top_tf_df = pd.read_csv(file_path, sep="\t")
                top_tfs = list(top_tf_df["TF"])
                tf_without_parentheses = [re.sub(r'\(.*?\)', '', s) for s in top_tfs]
                if tpm_filter > 0:
                    #top_tfs = [top_tf for top_tf in top_tfs if any(tpm_filtered_gen in top_tf for tpm_filtered_gen in tpm_filtered_valid_genes_hgnc)] # old  method
                    top_tfs = [top_tf for top_tf in tf_without_parentheses if all(tf_part in tpm_filtered_valid_genes_hgnc for tf_part in top_tf.split("::"))]
                if modus in ["list", "set" ]:                
                    if modus == "set":
                        top_tfs = set(top_tfs)
                    return top_tfs
                elif modus == "dataframe":
                    top_tf_df["TF"] = tf_without_parentheses
                    top_tf_df_filtered = top_tf_df[top_tf_df["TF"].isin(top_tfs)]
                    return top_tf_df_filtered


def create_overlap_df_with_tpm(comparison_main_dir, tpm_filter=0, modus: Union["list", "dataframe", "set"]="set"):
    """As input directory a main comparisonfodler should be given,
    e.g. /nfs/data3/IHEC/TF_PRIO/comparisons/cells. Based on the input folder
    it iterates through every subdirectory and if this subdirectory contains the top TF
    list files (dcg.tsv), it will read them out and, together with the celltypes and chromatin status,
    and return all TFs as lists as value and with the according name of the folder as key in a dictionary.
    """
    
    gtf_dict = read_gtf_and_create_ensg_hgnc_dict()
    dir_list = [path.path for path in os.scandir(comparison_main_dir) if path.is_dir()]
    overlap_df = pd.DataFrame(columns=['suppressive_marks','active_marks'])
    for dir_path in dir_list:
        comparison_name, chromatin_status = get_involved_celltypes_and_chromatin_status(dir_path)
        if chromatin_status == "custom-HM-combination":
            continue
        
        top_tf_comparison_name = "-".join([comparison_name, chromatin_status])
        top_tfs = read_out_top_tf_from_comparison_folder_with_tpm(dir_path, gtf_dict=gtf_dict, tpm_filter=tpm_filter, modus=modus)
        if top_tfs is not None:
            if not comparison_name in overlap_df.index.tolist():
                overlap_df.loc[comparison_name] = [None, None]
            if top_tfs is None:
                top_tfs = set()
            overlap_df.at[comparison_name, chromatin_status] = top_tfs

    return overlap_df


def read_out_top_tf_hm_from_comparison_folder_tpm_filter(comparison_path: Union[str, Path], gtf_dict, tpm_filter=0, modus: Union["dataframe", "set"]="set"):
    """Takes a path, checks it if it contains a directory with the top  TFs of the TF Prio run,
    and if so, reads the top TFs out and returns them as a list.
    """
    top_tf_hm_dict = {}
    #print(comparison_path)
    tpm_filtered_valid_genes_hgnc = []
    tpm_filtered_valid_genes_hgnc = read_out_tpm_filtered_genes(comparison_dir=comparison_path, gtf_dict=gtf_dict, tpm_filter_value=tpm_filter, as_hgnc=True, as_united_list=True)
    top_tf_hm_dir = os.path.join(comparison_path, "output/output/org_exbio_tfprio_steps_distributionAnalysis_CalculateDcgPerHm/output/")
    if not os.path.isdir(top_tf_hm_dir):
        print("The DcgRankPerHm directory "+str(top_tf_hm_dir)+" does not exist")
        return None
    else:
        if not os.listdir(top_tf_hm_dir):
            print("The DcgRankPerHm directory "+str(top_tf_hm_dir)+" is empty")
            return None
        else:
            for top_tf_hm_file in os.listdir(top_tf_hm_dir):
                hm_name = top_tf_hm_file.split(".")[0]
                top_tf_hm_file_path = os.path.join(top_tf_hm_dir, top_tf_hm_file)
                if os.stat(top_tf_hm_file_path).st_size > 0:
                    hm_tf_df = pd.read_csv(top_tf_hm_file_path, sep="\t", names = ["TF", "SCORE"], header=None)
                    top_hm_tfs = list(hm_tf_df["TF"])
                    tf_without_parentheses = [re.sub(r'\(.*?\)', '', s) for s in top_hm_tfs]
                    if tpm_filter > 0:
                        top_hm_tfs = [top_tf for top_tf in tf_without_parentheses if all(tf_part in tpm_filtered_valid_genes_hgnc for tf_part in top_tf.split("::"))]
                    if modus == "set":       
                        top_hm_tfs = set(top_hm_tfs)
                        top_tf_hm_dict[hm_name] = top_hm_tfs
                    else:
                        hm_tf_df["TF"] = tf_without_parentheses
                        hm_tf_df_filtered = hm_tf_df[hm_tf_df["TF"].isin(top_hm_tfs)]
                        top_tf_hm_dict[hm_name] = hm_tf_df_filtered

                else:
                    top_tf_hm_dict[hm_name] = set()
                
    return top_tf_hm_dict


def read_out_top_tfs_per_hm_tpm_filter(comparison_main_dir, tpm_filter=0, modus: Union["dataframe", "set"]="set"):
    """As input directory a main comparisonfodler should be given,
    e.g. /nfs/data3/IHEC/TF_PRIO/comparisons/cells. Based on the input folder
    it iterates through every subdirectory and if this subdirectory contains the top TF
    list files (dcg.tsv), it will read them out and, together with the celltypes and chromatin status,
    and return all TFs as lists as value and with the according name of the folder as key in a dictionary.
    """
    gtf_dict = read_gtf_and_create_ensg_hgnc_dict()
    dir_list = [path.path for path in os.scandir(comparison_main_dir) if path.is_dir()]
    tf_hm_df = pd.DataFrame(columns=['H3K4me1','H3K4me3','H3K27ac', 'H3K36me3', 'H3K27me3', 'H3K9me3'])
    for dir_path in dir_list:
        comparison_name, chromatin_status = get_involved_celltypes_and_chromatin_status(dir_path)
        if chromatin_status == "custom-HM-combination":
            continue

        top_tfs_dict = read_out_top_tf_hm_from_comparison_folder_tpm_filter(dir_path, gtf_dict = gtf_dict, tpm_filter = tpm_filter, modus = modus)
        if top_tfs_dict is not None:
            if not comparison_name in tf_hm_df.index.tolist():
                tf_hm_df.loc[comparison_name] = [None, None, None, None, None, None]
            if chromatin_status == "suppressive_marks":
                top_tfs_dict.pop("H3K4me3", None)
            for hm_name,  top_tfs in top_tfs_dict.items():
                if top_tfs is not None:
                    tf_hm_df.at[comparison_name, hm_name] = top_tfs
                else:
                    tf_hm_df.at[comparison_name, hm_name] = set()

    return tf_hm_df


def upset_plot(module_map, title_extension="",min_subset_size=0, max_subset_size=0):
    
    title = "UpSet"
    if title_extension:
        title = f"{title} - {title_extension}"
        
    tools, modules = zip(*module_map.items())
    all_elems = list(set().union(*modules))
    df = pd.DataFrame([[e in st for st in modules] for e in all_elems], columns = tools)
    
    upset_df = df.groupby(list(tools)).size()
    if max_subset_size == 0:
        usp.plot(upset_df, sort_by="cardinality", show_counts="%d", min_subset_size=min_subset_size)
    else:
        usp.plot(upset_df, sort_by="cardinality", show_counts="%d", max_subset_size=max_subset_size, min_subset_size=min_subset_size)
    plt.suptitle(title)
    plt.show()


def get_tpm_filtered_genes_df(all_comparisons_path: Union[str, Path], tpm_filter=1):
    """
    Takes the path of all comparisons as input and creates a dataframe which has for each celltype all tpm filtered genes.
    The filtered genes, originally saved as ENSG, are translated to HGNC symbols (using a created ensg-hgnc dictionary).
    Because for each comparison the same RNA-seq data is used for each celltype, its not necessary to save the TPM-filtered genes for each comparison.
    Default TPM filter value is set to 1.
    """
    # Creating gtf_dictionary for translating ENSG symbols into HGNC.
    gtf_dict = read_gtf_and_create_ensg_hgnc_dict()
    tpm_dict = {"celltype":[], "TPM_filtered_genes":[]}
    dir_list = [path.path for path in os.scandir(all_comparisons_path) if path.is_dir()]
    
    for dir_path in dir_list:
        comparison_name, chromatin_status = get_involved_celltypes_and_chromatin_status(dir_path)
        if chromatin_status == "custom-HM-combination":
            continue
        comparison_path = os.path.join(all_comparisons_path, dir_path)
        tpm_filtered_valid_genes_hgnc_dict = read_out_tpm_filtered_genes(comparison_dir=comparison_path, gtf_dict=gtf_dict, 
                                                                    tpm_filter_value=tpm_filter, as_hgnc=True, as_united_list=False)
        if tpm_filtered_valid_genes_hgnc_dict is not None:
            for celltype, tpm_genes in tpm_filtered_valid_genes_hgnc_dict.items():
                if celltype not in tpm_dict["celltype"]:
                    tpm_dict["celltype"].append(celltype)
                    tpm_dict["TPM_filtered_genes"].append(tpm_genes)
                
    tpm_df = pd.DataFrame(tpm_dict).set_index("celltype")
    
    return tpm_df


def get_expression_df(all_comparisons_path: Union[str, Path], modus: Union["mean", "single_samples"]):
    """
    Takes the path of all comparisons as input and creates a dataframe which has for each celltype the mean expression for each gene.
    The filtered genes, originally saved as ENSG, are translated to HGNC symbols (using a created ensg-hgnc dictionary).
    Because for each comparison the same RNA-seq data is used for each celltype, its not necessary to save the mean expression data for each comparison.
    """

    # Creating gtf_dictionary for translating ENSG symbols into HGNC.
    gtf_dict = read_gtf_and_create_ensg_hgnc_dict()
    # create empty dataframe to append columns with values to
    expression_df = pd.DataFrame({})
    # create list of all directories in the given input path
    dir_list = [path.path for path in os.scandir(all_comparisons_path) if path.is_dir()]
    
    for dir_path in dir_list:
        comparison_path = os.path.join(all_comparisons_path, dir_path)
        comparison_expression_dict = read_out_genes_expression(comparison_dir=comparison_path,  modus=modus)
        if comparison_expression_dict is not None:
            for celltype, celltype_expression_df in comparison_expression_dict.items():
                if modus == "mean":
                    if "gene" not in expression_df.columns:
                        expression_df["gene"]  =  celltype_expression_df["gene_id"]
                    if celltype not in expression_df.columns:
                        expression_df[celltype] = celltype_expression_df["mean"]  
                else: #single samples
                    if "gene" not in expression_df.columns:
                        expression_df["gene"]  =  celltype_expression_df["gene_id"]
                    if not any(celltype in col for col in expression_df.columns):
                        expression_df = pd.concat([expression_df, celltype_expression_df.iloc[:, 1:]], axis=1)

    
    #convert the ENSG symbols to HGNC
    expression_df["gene"] = [gtf_dict.get(symbol, None) for symbol in expression_df["gene"]]
    
    return expression_df


def read_out_genes_expression(comparison_dir, modus: Union["mean", "single_samples"]):
    """
    Reads each calculated mean expression file for a given comparison path in and saves the dataframe as value and, with the name of the mean expression file
    (therefor the celltype/ condition name) as key, in a dictionary and gives the dictionary back.
    """
    mean_expression_genes_dict = {}
    if modus == "mean":
        tpm_dir = os.path.join(comparison_dir,"output/output/org_exbio_tfprio_steps_rnaSeq_MeanExpression/output/")
    else:
        tpm_dir = os.path.join(comparison_dir,"output/output/org_exbio_tfprio_steps_rnaSeq_MeanExpression/input/")
    if not os.path.isdir(tpm_dir):
        print("The " +modus+ " expression directory "+str(tpm_dir)+" does not exist")
        return mean_expression_genes_dict
    else:
        if not os.listdir(tpm_dir):
            print("The " +modus+ " expression directory "+str(tpm_dir)+" is empty")
            return mean_expression_genes_dict
        else:
            for cell_type_file in os.listdir(tpm_dir):
                cell_type = cell_type_file.split(".")[0]
                file_path = os.path.join(tpm_dir, cell_type_file)
                file = pd.read_csv(file_path, sep="\t")
                if modus == "single_samples":
                    file.columns = [file.columns[0]] + [f"{cell_type}_{i}" for i in range(1,len(file.columns))]
                mean_expression_genes_dict[cell_type] = file
    
    return mean_expression_genes_dict


def create_top_tf_df_with_tpm(path_to_main_comparison_dir, tpm_filter=0):
    """Function for the creation of top TF df. Just with a more fitting name
    TODO: rename all occurences of the old function with the new one and after that change new one and delete old one
    """
    return create_overlap_df_with_tpm(comparison_main_dir, tpm_filter)


def read_out_generic_TFs(tf_df, occurence_limit = 9):
    """Takes a dataframe with the TFs as input. 
    The input df NEEDS to have 2 columns with "active_markers" and "repressive_markers" for the function to work.
    The input df should be created using the create_top_tf_df_with_tpm function.
    Returns two lists of TFs, that occur at least X times (set by the "occurence_limit" parameter) in all comparison
    """
    list_all_active_marks = [tf for tf_set in tf_df["active_marks"] for tf in tf_set]
    tf_counter_active_marks = Counter(list_all_active_marks)
    generic_TFs_active = [TF for TF, count in tf_counter_active_marks.items() if count >= occurence_limit ]

    #list_all_suppressive_markers = [tf for tf_set in tf_df["suppressive_markers"] for tf in tf_set]
    list_all_suppressive_marks = [tf for tf_set in tf_df["suppressive_marks"] for tf in tf_set]
    tf_counter_suppressive_marks = Counter(list_all_suppressive_marks)
    generic_TFs_suppressive = [TF for TF, count in tf_counter_suppressive_marks.items() if count >= occurence_limit ]

    return generic_TFs_active, generic_TFs_suppressive


def read_out_filtered_target_genes_per_celltype_dict(tf_df, all_comparisons_path, modus:Union["all_genes", "top_k_genes"]="top_k_genes", affinity_filter_value=0.05, tf_per_hm=False):
    """Reads out the target genes above the set affinity_filter_value. Either using the top_k_genes from the output folder or the files  from the input folder in 
    org_exbio_tfprio_steps_distributionAnalysis_TopKTargetGenes

    """
    tf_dict = tf_df.to_dict("index")
    tf_filtered_tg_dict={}
    dir_list = [path.path for path in os.scandir(all_comparisons_path) if path.is_dir()]
    for comparison, values in tf_dict.items():
        if comparison not in tf_filtered_tg_dict.keys():
            tf_filtered_tg_dict[comparison] = {}
            for col_name in tf_df.columns:
                tf_filtered_tg_dict[comparison][col_name] = list()

        for dir in dir_list:
            if comparison in dir:
                comparison_path = os.path.join(all_comparisons_path, dir)
                base_dir = os.path.join(comparison_path, "output/output/org_exbio_tfprio_steps_distributionAnalysis_TopKTargetGenes")
                if modus == "top_k_genes":
                    for mode, tf_list in values.items():
                        tf_files_list = [tf+".tsv" for tf in tf_list]
                        gene_base_dir = os.path.join(base_dir, "output", comparison)
                        for root, dirs, files in os.walk(gene_base_dir):
                            for tf_file in files:
                                if tf_file in tf_files_list:
                                    #print(tf_file)
                                    tf_file_path = os.path.join(root, tf_file)
                                    file = pd.read_csv(tf_file_path, sep="\t")
                                    file = file[file["Affinity"]>=affinity_filter_value]
                                    new_tgs = file["Gene"].tolist()
                                    tf_filtered_tg_dict[comparison][mode].extend(new_tgs)
            
                elif modus == "all_genes":
                    gene_base_dir = os.path.join(base_dir, "input", comparison)
                    for root, dirs, files in os.walk(gene_base_dir):
                        for tg_matrix in files:
                            if tg_matrix == "dcg.tsv":
                                continue
                            tg_matrix_path = os.path.join(root, tg_matrix)
                            tg_matrix_df = pd.read_csv(tg_matrix_path,  sep="\t")
                            for mode, tf_list in values.items():
                                for tf in tf_list:
                                    if tf in tg_matrix_df.columns:
                                        new_tgs = tg_matrix_df[tg_matrix_df[tf]>=affinity_filter_value]["Gene"].tolist()
                                        tf_filtered_tg_dict[comparison][mode].extend(new_tgs)
                                        # if consensus TGs:
                                        #tf_filtered_tg_dict[comparison][mode].intersection(new_tgs)
                            #for tf in tf_list:
                            #    if tf in tg_matrix_df.columns:
                            #        new_tgs = tg_matrix_df[tg_matrix_df[tf]>=affinity_filter_value]["Gene"].tolist()
                            #        list_of_filtered_genes.extend(new_tgs)

    return tf_filtered_tg_dict


def read_out_master_regulator_tf(tf_tgs_df, tf_df, all_comparisons_path, min_num_tg_total = 20000, min_num_tgs_per_tf = 1000, affinity_filter_value=0.05, per_cell_type=False):
    """Reads out the target genes above the set affinity_filter_value. Either using the top_k_genes df from the output folder or the files  from the input folder in 
    org_exbio_tfprio_steps_distributionAnalysis_TopKTargetGenes
    Needs as input 2 dataframes:
    tf_tgs_df should contain the dataframe with all target genes for the TFs of that specific condition
    tf_df should be the dataframe  which was used to create the tf_tgs_df
    "per_cell_type" set to True when the input dfs have not e.g. "macrophage-healthy-monocyte-healthy" but "macrophage-healthy"

    """
    tf_tg_dict = tf_tgs_df.to_dict("index")
    tf_filtered_tg_dict={}
    dir_list = [path.path for path in os.scandir(all_comparisons_path) if path.is_dir()]
    read_counts = 0
    for comparison, columns in tf_tg_dict.items():
        for dir in dir_list:
            if comparison in dir:
                comparison_path = os.path.join(all_comparisons_path, dir)
                if per_cell_type:
                    # comparison in this case is just the cell type
                    gene_base_dir = os.path.join(comparison_path, "output/output/org_exbio_tfprio_steps_distributionAnalysis_TopKTargetGenes/input", comparison)
                else:
                    gene_base_dir = os.path.join(comparison_path, "output/output/org_exbio_tfprio_steps_distributionAnalysis_TopKTargetGenes/input")
                for root, dirs, files in os.walk(gene_base_dir):
                    for tg_matrix in files:
                        if tg_matrix == "dcg.tsv":
                            continue
                        tg_matrix_path = os.path.join(root, tg_matrix)
                        tg_matrix_df = pd.read_csv(tg_matrix_path,  sep="\t")
                        read_counts += 1
                        for hm_markers, tg_list in columns.items():
                            if hm_markers.isin(["active_markers", "active_marks"])  and any(hm_mark in tg_matrix for hm_mark in ["H3K9me3", "H3K27me3"]):
                                #print("no TGs for: "+hm_markers+" in file: "+tg_matrix)
                                continue
                            if hm_markers.isin(["suppressive_markers", "suppressive_marks"])  and any(hm_mark in tg_matrix for hm_mark in ["H3K4me1", "H3K4me3", "H3K27ac", "H3K36me3"]):
                                #print("no TGs for: "+hm_markers+" in file: "+tg_matrix)
                                continue
                            if len(tg_list) >= min_num_tg_total:
                                tf_filtered_tg_dict[comparison]={}
                                tf_filtered_tg_dict[comparison][hm_markers]={}
                                tf_list = tf_df.loc[comparison, hm_markers]
                                for tf in tf_list:
                                    new_tgs = tg_matrix_df[tg_matrix_df[tf]>=affinity_filter_value]["Gene"].tolist()
                                    if len(new_tgs) >= min_num_tgs_per_tf:
                                        if tf not in tf_filtered_tg_dict[comparison][hm_markers].keys():
                                            tf_filtered_tg_dict[comparison][hm_markers][tf] = set()
                                        if tf in tg_matrix_df.columns:
                                            tf_filtered_tg_dict[comparison][hm_markers][tf].update(new_tgs)


    master_regulator_tf_df = flatten_master_regulator_dict(tf_filtered_tg_dict)
    print(read_counts)
    return master_regulator_tf_df


def flatten_master_regulator_dict(tf_dict):
    # List to store flattened data
    flattened_data = []

    # Recursive function to flatten the nested dictionary
    def flatten(d, level1=None, level2=None):
        for key, value in d.items():
            if isinstance(value, dict):
                if level1 is None:
                    flatten(value, key)
                elif level2 is None:
                    flatten(value, level1, key)
                else:
                    flatten(value, level1, level2, key)
            else:
                # Append the flattened data to the list
                flattened_data.append({
                    "comparison": level1,
                    "markers": level2,
                    "tf": key,
                    "target_strings": list(value)
                })

    # Start flattening from the top level
    flatten(tf_dict)

    # Convert the list of dictionaries to a DataFrame
    df = pd.DataFrame(flattened_data)
    return df


def enrichr_adjust_y_axis_title_font_sizes(enr_bg_results_df):
    """ Function for adjusting the font sizes for the dotplots genered by enrichr result from gseapy
    """
    results_p_val_filtered = enr_bg_results_df[enr_bg_results_df["Adjusted P-value"]<0.05]
    top_rows_per_set = results_p_val_filtered.groupby("Gene_set").head(5)
    number_of_enriched_sets_per_col = top_rows_per_set["Gene_set"].value_counts().sort_index()
    total_number_of_displayed_enriched_sets = number_of_enriched_sets_per_col.sum()
    if total_number_of_displayed_enriched_sets > 15:
        y_axis_size = 9
        title_size = 15
    else:
        y_axis_size = 12
        title_size = 20
    return y_axis_size, title_size


def enrichr_adjust_axis_title_font_sizes(enr_bg_results_df):
    """ Function for adjusting the font sizes for the dotplots genered by enrichr result from gseapy
    """
    results_p_val_filtered = enr_bg_results_df[enr_bg_results_df["Adjusted P-value"]<0.05]
    top_rows_per_set = results_p_val_filtered.groupby("Gene_set").head(5)
    number_of_enriched_sets_per_col = top_rows_per_set["Gene_set"].value_counts().sort_index()
    total_number_of_displayed_enriched_sets = number_of_enriched_sets_per_col.sum()
    if total_number_of_displayed_enriched_sets > 25:
        y_axis_size = 8
        x_axis_size = 9
        title_size = 15
    elif total_number_of_displayed_enriched_sets > 15:
        y_axis_size = 9
        x_axis_size = 9
        title_size = 15
    else:
        y_axis_size = 12
        x_axis_size = 12
        title_size = 20
    return x_axis_size, y_axis_size, title_size


def read_in_gene_ontology_file(file):
    """Reads in GO term files. Adjusts for the shift in the column names,
    e.g. first column name is "index" and the extra NaN value column as the last columns as a result of
    the col name shift. Also changes the the HGNC symbol to upper to make it easier comparable to other genes.
    """
    go_info = pd.read_csv(file, sep="\t").reset_index()
    go_info.columns.values[:11] = go_info.columns.values[1:]
    go_info = go_info.iloc[:, :-1]
    go_info["Symbol"] = go_info["Symbol"].str.upper()
    return go_info
