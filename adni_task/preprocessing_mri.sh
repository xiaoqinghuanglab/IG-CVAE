#!/bin/bash

# Check if an input directory is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input_directory>"
    exit 1
fi

INPUT_DIR=$1
OUTPUT_DIR="MRI_Preprocessed_Data"
INTERMEDIATE_DIR="${OUTPUT_DIR}/intermediate"

# Create output and intermediate directories if they don't exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$INTERMEDIATE_DIR"

# Loop through all NIfTI images in the input directory
for img in "${INPUT_DIR}"/*.nii; do
    if [[ ! -f $img ]]; then
        echo "No NIfTI files found in ${INPUT_DIR}."
        exit 1
    fi

    filename=$(basename -- "$img")
    base="${filename%.nii}"

    echo "Processing $filename..."

    # Intermediate file paths
    tmp_skullstripped="${INTERMEDIATE_DIR}/${base}_skullstripped.nii"
    tmp_reoriented="${INTERMEDIATE_DIR}/${base}_reoriented.nii"
    tmp_registered="${INTERMEDIATE_DIR}/${base}_registered.nii"
    tmp_biascorrected="${INTERMEDIATE_DIR}/${base}_biascorrected.nii"
    tmp_matrix="${INTERMEDIATE_DIR}/${base}_MRI_to_MNI.mat"

    # Step 1: Skull stripping
    bet2 "$img" "$tmp_skullstripped" -f 0.5

    # Step 2: Reorient to standard orientation
    fslreorient2std "$tmp_skullstripped" "$tmp_reoriented"

    # Step 3: Register to MNI152 space
    flirt -in "$tmp_reoriented" -ref $FSLDIR/data/standard/MNI152_T1_2mm.nii.gz -dof 12 -omat "$tmp_matrix" -out "$tmp_registered"

    # Step 4: Bias field correction
    fast -B "$tmp_registered"
    mv "${tmp_registered%.nii}_restore.nii.gz" "${tmp_biascorrected}.gz"
    gunzip -f "${tmp_biascorrected}.gz"

    # Step 5: Intensity normalization
    mean=$(fslstats "$tmp_biascorrected" -M)
    std=$(fslstats "$tmp_biascorrected" -S)
    final_output="${OUTPUT_DIR}/${base}.nii"
    fslmaths "$tmp_biascorrected" -sub $mean -div $std "$final_output"

    echo "Finished processing $filename. Final output: ${final_output}"
done

echo "✅ Batch processing complete! Intermediate files saved to: ${INTERMEDIATE_DIR}"

