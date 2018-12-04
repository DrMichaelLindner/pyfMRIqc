# fMRI_QC
fMRI_QC.py calculates and provides information of a given functional MRI nifti file for a quality check.


## *Usage*
        python fmri_qc.py 
        python fmri_qc.py -n <func_nift_file> 
        python fmri_qc.py -n <func_nift_file> -m <motion_file> 
        python fmri_qc.py -n <func_nift_file> -m <motion_file> -t <mask_threshold>
        python fmri_qc.py -n <func_nift_file> -m <motion_file> -t <mask_threshold> -o <output_folder>
        
        
## *Input*
        -n:   functional MR nifti file 
        -m:   motion parameters file of motion correction from FSL (*.par) or SPM (rp*.txt)
        -t:   threshold of mean values for the mask to calculate SNR etc.
        -o:   output folder  
    All input are optionally: if not defined input dialogs will pop up to select the files manually.
    
    
## *Output*
fMRI_QC.py creates the following outputs:

nifti files
    - MEAN over time
    - VAR over time
    - MASK (binary - containing voxels above the threshold)
    - SNR signal-to-noise ratio
    - SQUARED DIFF
    - SQUARED SCALED DIFF (Squared Diff / Mean Diff)

a png image containing the following plots:
    - scaled variability: Mean (over all voxel) squared diff plot over time / global mean squared diff
    - slice by slice variability: Mean (mean per slice) squared diff plot over time / global mean squared diff
    - if motion file was selected: sum of relative movements over time (z-scored)
    - scaled mean voxel intensity: mean(data/global mean) )(z-scored)
    - variance of scaled variability (z-scored)
    - min/mean/max slice variability

a text file
    -  containing an overview about scan and QC parameters

See guidelines about the output in the [Wiki](https://github.com/DrMichaelLindner/fMRI_QC/wiki)


## *Install*  
Copy the fMRI_QC folder in a folder of your choice on your system and add the directory to your PYTHONPATH.


## *Dependencies*  
fMRI_QC is developed in python 3.6 and the following packages need to be installed: 
numpy, nibabel, matplotloib, scipy and easygui. 
(fMRI_QC was developed and tested with the following versions: 
python: 3.6.4, numpy: 1.14.3, nibabel: 2.2.1, matplotloib: 2.1.0, scipy: 0.19.1 and easygui: 0.98.1 )

    
## *License*  
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License (GPLv3) as published
by the Free Software Foundation;

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
  
  
## *Author*
Michael Lindner  
University of Reading, 2018  
School of Psychology and Clinical Language Sciences  
Centre for Integrative Neuroscience and Neurodynamics
