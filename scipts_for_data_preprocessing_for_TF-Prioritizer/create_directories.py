#!/usr/bin/env python3
import os
import numpy
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--input_dir', required=True)
parser.add_argument('--biomaterial_type', choices=['tissue', 'cells'], required=True)

parser.add_argument("--sample1", choices=["blood-leukemia", "brain-healthy", "brain-autism",
                                          "t-cell-healthy", "neutrophil-healthy", "monocyte-healthy", "b-cells-healthy",
                                          "b-cells-leukemia", "myeloid-cells-leukemia", "macrophage-healthy"])
parser.add_argument("--sample2", choices=["blood-leukemia", "brain-healthy", "brain-autism",
                                          "t-cell-healthy", "neutrophil-healthy", "monocyte-healthy", "b-cells-healthy",
                                          "b-cells-leukemia", "myeloid-cells-leukemia", "macrophage-healthy"])

hist_mods = parser.add_mutually_exclusive_group(required=True)
hist_mods.add_argument('--histone_modifications_single_choice', nargs="+", choices=["h3k27me3", "H3K9me3", "H3K4me3","H3K4me1", "H3k27ac", "h3k36me3"])
hist_mods.add_argument('--histone_modifications_modus', choices=["heterochromatin", "euchromatin"])
parser.add_argument('--mrna', action='store_true')
parser.add_argument('--total_rna', action='store_true')

args = parser.parse_args()

tissue_samples = ["blood-leukemia", "brain-healthy", "brain-autism"]
cells_samples = ["t-cell-healthy", "neutrophil-healthy", "monocyte-healthy", "b-cells-healthy",
                      "b-cells-leukemia", "myeloid-cells-leukemia", "macrophage-healthy"]

basedir = args.input_dir
samples = [args.sample1, args.sample2]
biomaterial = args.biomaterial_type
h3k27me3_val = False
h3k27ac_val = False
h3k9me3_val = False
h3k4me3_val = False
h3K4me1_val = False
h3k36me3_val = False
mrna_val = False
total_rna_val = False

#print("mrna:" + str(mrna_val))
#print("total_rna:" + str(total_rna_val))
#print(biomaterial)

if args.mrna is not None:
    mrna_val = args.mrna
if args.total_rna is not None:
    total_rna_val = args.total_rna
#print("mrna:" + str(mrna_val))
#print("total_rna:" + str(total_rna_val))

#if args.biomaterial_type == "tissue": #biomaterial
if biomaterial == "tissue":
    if args.sample1 in cells_samples or args.sample2 in cells_samples:
        parser.error("Both samples have to be in the chosen biomaterial type")
if biomaterial == "cells":
    if args.sample1 in tissue_samples or args.sample2 in tissue_samples:
        parser.error("Both samples have to be in the chosen biomaterial type")

if args.sample1 == args.sample2:
    parser.error("sample1 and sample2 have to be different samples")

if args.sample1 in tissue_samples and args.sample2 in cells_samples or args.sample1 in cells_samples and args.sample2 in tissue_samples:
    parser.error("sample1 and sample2 have to be the same biomaterial type")


if args.histone_modifications_modus is not None:
    if args.histone_modifications_modus == "heterochromatin":
        h3k27me3_val = True
        h3k9me3_val = True
        h3k4me3_val = True
    if args.histone_modifications_modus == "euchromatin":
        h3k4me3_val = True
        h3K4me1_val = True
        h3k27ac_val = True
        h3k36me3_val = True
if args.histone_modifications_single_choice is not None:
    if "h3k27me3" in args.histone_modifications_single_choice:
        h3k27me3_val = True
    if "h3k27ac" in args.histone_modifications_single_choice:
        h3k27ac_val = True
    if "h3k9me3" in args.histone_modifications_single_choice:
        h3k9me3_val = True
    if "h3k4me3" in args.histone_modifications_single_choice:
        h3k4me3_val = True
    if "h3K4me1" in args.histone_modifications_single_choice:
        h3K4me1_val = True
    if "h3k36me3" in args.histone_modifications_single_choice:
        h3k36me3_val = True


def create_dir_and_symlinks(basedir, biomaterial_type, samplelist, h3k27ac, h3k27me3, h3k36me3, h3k4me1, h3k4me3, h3k9me3, m_rna, total_rna):
    basedir = "/nfs/data3/IHEC/TF_PRIO/comparisons"
    main_dir_path = create_main_dir_path(basedir = basedir, samplelist = samplelist, biomaterial_type=biomaterial_type, h3k27ac = h3k27ac, h3k27me3 = h3k27me3, h3k36me3 = h3k36me3,
                                     h3k4me1 = h3k4me1, h3k4me3 = h3k4me3, h3k9me3 = h3k9me3, m_rna = m_rna, total_rna = total_rna)
    create_directories(main_dir_path, samplelist)
    create_symlinks(main_dir_path, samplelist, h3k27ac = h3k27ac, h3k27me3 = h3k27me3, h3k36me3 = h3k36me3,
                                     h3k4me1 = h3k4me1, h3k4me3 = h3k4me3, h3k9me3 = h3k9me3, mrna = m_rna, total_rna = total_rna)
    delete_empty_directories(main_dir_path)


def create_main_dir_path(basedir, samplelist, biomaterial_type = "", h3k27ac = False, h3k27me3 = False, h3k36me3 = False,
                     h3k4me1 = False, h3k4me3 = False, h3k9me3 = False, m_rna = False, total_rna = False):
    basedir = basedir
    if biomaterial_type == "tissue":
        basedir = os.path.join(basedir, "tissue")
    if biomaterial_type == "cells":
        basedir = os.path.join(basedir, "cells")

    dirname = ""
    for sample in samplelist:
        if dirname:
            dirname = dirname + "-" + sample
        else:
            dirname = sample
    if h3k27ac:
        dirname = dirname + "-H3K27ac"
    if h3k27me3:
        dirname = dirname + "-H3K27me3"
    if h3k36me3:
        dirname = dirname + "-H3k36me3"
    if h3k4me1:
        dirname = dirname + "-H3K4me1"
    if h3k4me3:
        dirname = dirname + "-H3K4me3"
    if h3k9me3:
        dirname = dirname + "-H3K9me3"
    if m_rna:
        dirname = dirname + "-mRNA"
    if total_rna:
        dirname = dirname + "-totalRNA"
    main_dir = os.path.join(basedir, dirname)
    return main_dir


def create_directories(maindir, samplelist):
    chip_seq_dir = os.path.join(maindir, "chipseq")
    histon_mods = ["H3K36me3", "H3K27ac", "H3K27me3", "H3K4me1", "H3K4me3", "H3K9me3"]
    if not os.path.exists(chip_seq_dir):
        os.makedirs(chip_seq_dir)

    for sample in samplelist:
        chip_seq_dir_sample = os.path.join(chip_seq_dir, sample)
        if not os.path.exists(chip_seq_dir_sample):
            os.mkdir(chip_seq_dir_sample)
        for histon_mod in histon_mods:
            chip_seq_dir_sample_histon = os.path.join(chip_seq_dir_sample, histon_mod)
            os.mkdir(chip_seq_dir_sample_histon)

    rna_seq_dir = os.path.join(maindir, "rna-seq")
    if not os.path.exists(rna_seq_dir):
        os.makedirs(rna_seq_dir)

    for sample in samplelist:
        rna_seq_dir_sample = os.path.join(rna_seq_dir, sample)
        if not os.path.exists(rna_seq_dir_sample):
            os.mkdir(rna_seq_dir_sample)


def create_symlinks(main_dir_path, samplelist, h3k27ac, h3k27me3, h3k36me3, h3k4me1, h3k4me3, h3k9me3, mrna, total_rna):
    ### uses the df_path and epirr_id functions and gives the values
    path_df = create_path_dataframe()
    epirr_id_dict = get_precast_epirr_ids(path_df)
    for sample in samplelist:
        sample_epirr_ids = epirr_id_dict[sample]
        write_chip_seq_symlinks(sample_epirr_ids, path_df, os.path.join(main_dir_path, "chipseq", sample),
                                h3k27ac = h3k27ac, h3k27me3 = h3k27me3, h3k36me3 = h3k36me3, h3k4me1 = h3k4me1, h3k4me3 = h3k4me3, h3k9me3 = h3k9me3)
        write_rna_seq_symlinks(sample_epirr_ids, path_df, os.path.join(main_dir_path, "rna-seq", sample), mrna = mrna, total_rna = total_rna)


#def create_path_dataframe(path_to_epiatlas, path_to_IHEc_metadata):
def create_path_dataframe():
    # create new path df for paths of input files on server
    df_epiatlas = pd.read_csv("epiatlas_metadata.csv")
    df_metadata = pd.read_csv("IHEC_metadata_harmonization.v1.1.csv")
    df_metadata = df_metadata.replace("primary cell culture","primary cell")

    # filter metadata for sample_ontology and biomaterial_type
    df_metadata_filtered = df_metadata[["epirr_id_without_version","harmonized_biomaterial_type","harmonized_sample_ontology_intermediate","harmonized_sample_disease_high", "harmonized_sample_disease_intermediate","harmonized_sample_disease"]]

    # merge dataframes
    df_merge = pd.merge(df_epiatlas, df_metadata_filtered, left_on="epirr_id_without_version", right_on="epirr_id_without_version")

    # Filter merge df for histon modifications and just keep informative columns
    experiment_types_2_keep = ["H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K9me3","total-RNA-Seq","mRNA-Seq"]
    df_merge_filter = df_merge[df_merge["experiment_type"].isin(experiment_types_2_keep)]
    df_merge_filter = df_merge_filter[["epirr_id","epirr_id_without_version", "experiment_type", "data_file_path",
                                   "harmonized_biomaterial_type","harmonized_sample_ontology_intermediate",
                                   "harmonized_sample_disease_high","harmonized_sample_disease_intermediate","harmonized_sample_disease"]]

    # Cast long to wide
    df_path_links = df_merge_filter.pivot(index=["epirr_id","epirr_id_without_version","harmonized_sample_ontology_intermediate"
                                                 ,"harmonized_biomaterial_type","harmonized_sample_disease_high","harmonized_sample_disease_intermediate","harmonized_sample_disease"],
                                          columns="experiment_type", values="data_file_path")

    # Read in all available imputed peaks and fill in accordingly the paths if NaN values in path_df
    imputed_files_list = []
    with open("all_files.txt", "r") as file:
        for line in file:
            # Remove the trailing newline character and add the line to the list
            imputed_files_list.append(line.strip())

    for i in range(len(imputed_files_list)):
        new_val = imputed_files_list[i].split("/")[-1]
        imputed_files_list[i] = new_val

    for row in range(df_path_links.shape[0]):
        for col in range(df_path_links.shape[1]):
            if pd.isna(df_path_links.iat[row, col]):
                hist_mod = df_path_links.columns[col]
                epirr_id = df_path_links.index[row][0]
                file_name_list = "impute_"+epirr_id+"_"+hist_mod+".pval.bw.narrowPeak.gz"
                if file_name_list in imputed_files_list:
                    file_name = "impute_"+epirr_id+"_"+hist_mod+".pval.bw.narrowPeak"
                    path = "IHEC/incoming/ChIP-Seq_imputed/imputed_peaks/narrowPeak/"+hist_mod+"/"+file_name
                    df_path_links.iat[row, col] = path

    def change_path_suffix(value, suffix, new_suffix=".pval0.01.500K.narrowPeak"):
        if isinstance(value, str):
            if value.endswith(suffix):
                return value.replace(suffix, new_suffix)
            else:
                return value

    def change_path(value, new_path):
        if isinstance(value, str):
            file = value.split("/")[-1]
            new_file_path = os.path.join(new_path, file)
            return new_file_path

    df_index = df_path_links.reset_index()

    chip_seq_data = ["H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K9me3"]
    for chip_seq_experiment_type in chip_seq_data:
        df_index[chip_seq_experiment_type] = df_index[chip_seq_experiment_type].apply(lambda x: change_path_suffix(x, ".*", new_suffix=".pval0.01.500K.narrowPeak"))
        df_index[chip_seq_experiment_type] = df_index[chip_seq_experiment_type].apply(lambda x: change_path(x, new_path = "/nfs/data3/IHEC/ChIP-Seq_peaks"))

    rna_seq_data =["total-RNA-Seq","mRNA-Seq"]
    for rna_seq_experiment_type in rna_seq_data:
        df_index[rna_seq_experiment_type] = df_index[rna_seq_experiment_type].apply(lambda x: change_path_suffix(x, ".*", new_suffix=".tsv"))
        df_index[rna_seq_experiment_type] = df_index[rna_seq_experiment_type].apply(lambda x: change_path(x, new_path = "/nfs/data3/IHEC/TF_PRIO/RNA_files"))

    # remove all rows that are just mRNA and  total-RNA seq samples
    df_index = df_index.dropna(subset=["H3K27ac"])

    return df_index


def get_precast_epirr_ids(df_paths):
    epirr_id_dict = {}

    epirr_id_dict["blood-leukemia"] = get_epirr_ids(df_paths, "primary tissue", harmonized_sample_ontology_intermediate="venous blood",
                                                             harmonized_sample_disease_intermediate = "Leukemia", harmonized_sample_disease = "Chronic Lymphocytic Leukemia")
    epirr_id_dict["brain-healthy"] = get_epirr_ids(df_paths, "primary tissue", harmonized_sample_ontology_intermediate="brain", harmonized_sample_disease_high="Healthy/None")
    epirr_id_dict["brain-autism"] = get_epirr_ids(df_paths, "primary tissue", harmonized_sample_ontology_intermediate="brain",
                                        harmonized_sample_disease_high="Disease", harmonized_sample_disease_intermediate = "Psychiatric Disorder")
    epirr_id_dict["t-cell-healthy"] = get_epirr_ids(df_paths, "primary cell", harmonized_sample_ontology_intermediate="T cell", harmonized_sample_disease_high="Healthy/None")
    epirr_id_dict["neutrophil-healthy"] = get_epirr_ids(df_paths, "primary cell", harmonized_sample_ontology_intermediate="neutrophil", harmonized_sample_disease_high="Healthy/None")
    epirr_id_dict["monocyte-healthy"] = get_epirr_ids(df_paths, "primary cell", harmonized_sample_ontology_intermediate="monocyte", harmonized_sample_disease_high="Healthy/None")
    epirr_id_dict["b-cells-healthy"] = get_epirr_ids(df_paths, "primary cell", harmonized_sample_ontology_intermediate="lymphocyte of B lineage", harmonized_sample_disease_high="Healthy/None")
    epirr_id_dict["b-cells-leukemia"] = get_epirr_ids(df_paths, "primary cell", harmonized_sample_ontology_intermediate="lymphocyte of B lineage", harmonized_sample_disease_intermediate = "Leukemia")
    epirr_id_dict["myeloid-cells-leukemia"] = get_epirr_ids(df_paths, "primary cell", harmonized_sample_ontology_intermediate="myeloid cell", harmonized_sample_disease_intermediate = "Leukemia")
    epirr_id_dict["macrophage-healthy"] = get_epirr_ids(df_paths, "primary cell", harmonized_sample_ontology_intermediate="macrophage", harmonized_sample_disease_high="Healthy/None")

    return epirr_id_dict


def get_epirr_ids(df_paths, biomaterial_type, harmonized_sample_ontology_intermediate = None, harmonized_sample_disease_high = None, harmonized_sample_disease_intermediate = None, harmonized_sample_disease = None):
    df_filtered = df_paths[(df_paths["harmonized_biomaterial_type"]==biomaterial_type)]

    if harmonized_sample_ontology_intermediate is not None:
        df_filtered = df_filtered[(df_filtered["harmonized_sample_ontology_intermediate"]==harmonized_sample_ontology_intermediate)]
    if harmonized_sample_disease_high is not None:
        df_filtered = df_filtered[(df_filtered["harmonized_sample_disease_high"]==harmonized_sample_disease_high)]
    if harmonized_sample_disease_intermediate is not None:
        df_filtered = df_filtered[(df_filtered["harmonized_sample_disease_intermediate"]==harmonized_sample_disease_intermediate)]
    if harmonized_sample_disease is not None:
        df_filtered = df_filtered[(df_filtered["harmonized_sample_disease"]==harmonized_sample_disease)]

    epirr_id_list = df_filtered["epirr_id_without_version"]

    return epirr_id_list


def write_chip_seq_symlinks(epirr_id_list, df_paths, outputdir, h3k27ac = False, h3k27me3 = False, h3k36me3 = False, h3k4me1 = False, h3k4me3 = False, h3k9me3 = False):
    for index, row in df_paths[df_paths["epirr_id_without_version"].isin(epirr_id_list)].iterrows():
        if h3k27ac:
            os.link(row["H3K27ac"], os.path.join(outputdir, "H3K27ac", row["H3K27ac"].split("/")[-1]))
        if h3k27me3:
            os.link(row["H3K27me3"], os.path.join(outputdir, "H3K27me3", row["H3K27me3"].split("/")[-1]))
        if h3k36me3:
            os.link(row["H3K36me3"], os.path.join(outputdir, "H3K36me3", row["H3K36me3"].split("/")[-1]))
        if h3k4me1:
            os.link(row["H3K4me1"], os.path.join(outputdir, "H3K4me1", row["H3K4me1"].split("/")[-1]))
        if h3k4me3:
            os.link(row["H3K4me3"], os.path.join(outputdir, "H3K4me3", row["H3K4me3"].split("/")[-1]))
        if h3k9me3:
            os.link(row["H3K9me3"], os.path.join(outputdir, "H3K9me3", row["H3K9me3"].split("/")[-1]))


def write_rna_seq_symlinks(epirr_id_list, df_paths, outputdir, mrna = False, total_rna = False):
    for index, row in df_paths[df_paths["epirr_id_without_version"].isin(epirr_id_list)].iterrows():
        if mrna and total_rna:
            if row["mRNA-Seq"] and row["total-RNA-Seq"]:
                os.link(row["mRNA-Seq"], os.path.join(outputdir, row["mRNA-Seq"].split("/")[-1]))
            elif row["mRNA-Seq"]:
                os.link(row["mRNA-Seq"], os.path.join(outputdir, row["mRNA-Seq"].split("/")[-1]))
            elif row["total-RNA-Seq"]:
                os.link(row["total-RNA-Seq"], os.path.join(outputdir, row["total-RNA-Seq"].split("/")[-1]))
        elif mrna:
            if row["mRNA-Seq"]:
                os.link(row["mRNA-Seq"], os.path.join(outputdir, row["mRNA-Seq"].split("/")[-1]))
        elif total_rna:
            if row["total-RNA-Seq"]:
                os.link(row["total-RNA-Seq"], os.path.join(outputdir, row["total-RNA-Seq"].split("/")[-1]))


def delete_empty_directories(base_dir):
    for root, dirs, files in os.walk(base_dir, topdown=False):
        for directory in dirs:
            dir_path = os.path.join(root, directory)
            
            if not os.listdir(dir_path):
                try:
                    os.rmdir(dir_path)
                    print(f"Deleted empty directory: {dir_path}")
                except OSError as e:
                    print(f"Error: {e}")
            

create_dir_and_symlinks(basedir=basedir, biomaterial_type=biomaterial, samplelist=samples, h3k27ac=h3k27ac_val, h3k27me3=h3k27me3_val,
                        h3k36me3=h3k36me3_val, h3k4me1=h3K4me1_val, h3k4me3=h3k4me3_val, h3k9me3=h3k9me3_val, m_rna=mrna_val, total_rna=total_rna_val)
