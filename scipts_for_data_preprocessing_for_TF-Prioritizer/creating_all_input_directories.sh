#!/bin/bash


biomaterial_options=("tissue" "cells")
sample_options=("blood-leukemia" "brain-healthy" "t-cell-healthy" "neutrophil-healthy" "monocyte-healthy" "b-cells-healthy" "macrophage-healthy" "b-cells-leukemia" "myeloid-cells""myeloid-cells-leukemia")
histon_modus_options=("heterochromatin" "euchromatin")

tissue_choices=${sample_options[@]:0:2}
cells_healthy=${sample_options[@]:2:5}
cells_leukemia=${sample_options[@]:7:2}

# Function to check if arguments 2 and 3 are different
check_samples() {
    if [ "$1" != "$2" ]; then
        return 0  # Arguments are different
    else
        return 1  # Arguments are the same
    fi
}

# cells - healthy - healthy inter comparison
for sample1 in $cells_healthy; do
    for sample2 in $cells_healthy; do
        if check_samples "$sample1" "$sample2"; then
            for histon_modus in "${histon_modus_options[@]}"; do
                python create_directories.py --input_dir "/nfs/proj/ihec_tfprio/comparisons" --biomaterial_type "cells" --sample1 "$sample1" --sample2 "$sample2" --histone_modifications_modus "$histon_modus" --mrna --total_rna
            done
        fi
    done
done

# cells - healthy - leukemia inter comparison
for sample1 in $cells_healthy; do
    for sample2 in $cells_leukemia; do
        if check_samples "$sample1" "$sample2"; then
            for histon_modus in "${histon_modus_options[@]}"; do
                python create_directories.py --input_dir "/nfs/proj/ihec_tfprio/comparisons" --biomaterial_type "cells" --sample1 "$sample1" --sample2 "$sample2" --histone_modifications_modus "$histon_modus" --mrna --total_rna
            done
        fi
    done
done

# cells - leukemia - leukemia inter comparison
for sample1 in $cells_leukemia; do
    for sample2 in $cells_leukemia; do
        if check_samples "$sample1" "$sample2"; then
            for histon_modus in "${histon_modus_options[@]}"; do
                python create_directories.py --input_dir "/nfs/proj/ihec_tfprio/comparisons" --biomaterial_type "cells" --sample1 "$sample1" --sample2 "$sample2" --histone_modifications_modus "$histon_modus" --mrna --total_rna
            done
        fi
    done
done

# cells - healthy - leukemia intra comparison
for histon_modus in "${histon_modus_options[@]}"; do
                python create_directories.py --input_dir "/nfs/proj/ihec_tfprio/comparisons" --biomaterial_type "cells" --sample1 "b-cells-healthy" --sample2 "b-cells-leukemia" --histone_modifications_modus "$histon_modus" --mrna --total_rna
done

#tissue - healthy - leukemia inter
for histon_modus in "${histon_modus_options[@]}"; do
                python create_directories.py --input_dir "/nfs/proj/ihec_tfprio/comparisons" --biomaterial_type "tissue" --sample1 "blood-leukemia" --sample2 "brain-healthy" --histone_modifications_modus "$histon_modus" --mrna --total_rna
done