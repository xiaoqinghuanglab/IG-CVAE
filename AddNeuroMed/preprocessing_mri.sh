#!/bin/bash

# Check if an input directory is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input_directory>"
    exit 1
fi

INPUT_DIR=$1
OUTPUT_DIR="MRI_Processed"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all NIfTI images in the input directory
for img in "${INPUT_DIR}"/*.nii.gz; do
    if [[ ! -f $img ]]; then
        echo "No NIfTI files found in ${INPUT_DIR}."
        exit 1
    fi

    filename=$(basename -- "$img")
    base="${filename%.nii.gz}"

    echo "Processing $filename..."

    # Temporary file names (saved in /tmp to avoid clutter)
    tmp_skullstripped="/tmp/${base}_skullstripped.nii.gz"
    tmp_reoriented="/tmp/${base}_reoriented.nii.gz"
    tmp_registered="/tmp/${base}_registered.nii.gz"
    tmp_biascorrected="/tmp/${base}_biascorrected.nii.gz"

    # Step 1: Skull stripping
    bet2 "$img" "$tmp_skullstripped" -f 0.5

    # Step 2: Reorient to standard orientation
    fslreorient2std "$tmp_skullstripped" "$tmp_reoriented"

    # Step 3: Register to MNI152 space
    flirt -in "$tmp_reoriented" -ref $FSLDIR/data/standard/MNI152_T1_2mm.nii.gz -dof 12 -omat /tmp/${base}_MRI_to_MNI.mat -out "$tmp_registered"

    # Step 4: Bias field correction
    fast -B "$tmp_registered"
    mv "${tmp_registered%.nii.gz}_restore.nii.gz" "$tmp_biascorrected"

    # Step 5: Intensity normalization
    mean=$(fslstats "$tmp_biascorrected" -M)
    std=$(fslstats "$tmp_biascorrected" -S)
    final_output_gz="${OUTPUT_DIR}/${base}.nii.gz"
    fslmaths "$tmp_biascorrected" -sub $mean -div $std "$final_output_gz"

    # Step 6: Convert .nii.gz to .nii
    gunzip -f "$final_output_gz"   # gunzip will replace .nii.gz → .nii

    echo "Finished processing $filename. Final output: ${OUTPUT_DIR}/${base}.nii"
    
    # Clean up temporary files
    rm -f "$tmp_skullstripped" "$tmp_reoriented" "$tmp_registered" "$tmp_biascorrected" /tmp/${base}_MRI_to_MNI.mat
done

echo "Batch processing complete!"

