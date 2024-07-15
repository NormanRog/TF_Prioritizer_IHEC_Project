import pandas as pd
import os
import json
import numpy as np

    
with open("/nfs/data3/IHEC/TF_PRIO/scripts/example_configs/chipSeq.json", "r") as file:
    config_chip = json.load(file)
    file.close()

base_dir = "/nfs/data3/IHEC/TF_PRIO/comparisons"
info_dir = "/nfs/data3/IHEC/TF_PRIO/additional_files"
chipAtlas_tissues = pd.read_csv("/nfs/data3/IHEC/TF_PRIO/additional_files/ChipAtlas_tissues.csv")

for biomaterial_type_dir_name in ["cells", "tissue"]:
    biomaterial_type_dir = os.path.join(base_dir, biomaterial_type_dir_name)
    for comparison_dir_name in [dir_name for dir_name in os.listdir(biomaterial_type_dir)]:
        comparison_main_dir = os.path.join(biomaterial_type_dir, comparison_dir_name)

        states_dict = {}
        tissues = []
        chip_seq_dir = os.path.join(comparison_main_dir,"chipseq","")
        for state in os.listdir(chip_seq_dir):
            states_dict[state] = state
            tissue = chipAtlas_tissues[chipAtlas_tissues["biomaterial_type"]==state]["chip_atlas_tissue"].values[0]
            tissues.append(tissue)

        tissues = list(set(tissues))
        config_chip["InputConfigs"]["peaks"] = chip_seq_dir
        config_chip["InputConfigs"]["rnaSeq"] = os.path.join(comparison_main_dir,"rna-seq","")
        config_chip["InputConfigs"]["geneIDs"] = os.path.join(info_dir, "Gene_id.txt")
        config_chip["InputConfigs"]["sameStages"] = states_dict
        config_chip["InputConfigs"]["geneAnnotationFile"] = os.path.join(info_dir, "gencode.v29.annotation.gtf")
        config_chip["InputConfigs"]["genome"] = "hg38"
        config_chip["InputConfigs"]["biomartSpecies"] = "hsapiens_gene_ensembl"
        config_chip["MixOptions"]["blackListPath"] = os.path.join(info_dir, "hg38-blacklist.v2.bed")
        config_chip["TEPIC"]["PWMs"] = os.path.join(info_dir, "combined_Jaspar_Hocomoco_Kellis_human_PSEM.PSEM")
        config_chip["TEPIC"]["referenceGenome"] = os.path.join(info_dir, "hg38_genome.fa")
        config_chip["TEPIC"]["onlyDNasePeaks"] = os.path.join(info_dir, "gencode.v29.annotation.gtf")
        config_chip["ChipAtlas"]["tissueTypes"] = tissues

        new_config_name = os.path.join(comparison_main_dir, "chipSeq.json")
        with open(new_config_name, 'w') as f:
            json.dump(config_chip, f, indent=2)
            f.close()
