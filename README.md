![LOGO](pyfMRIqc_logo.png)

# pyfMRIqc
pyfMRIqc.py calculates and provides information of a given functional MRI nifti file for a quality check.
pyfMRIqc.py can be used by giving input parameter or without input parameter. If no input parameter are defined
pyfMRIqc.py will guide the user through input dialogs to manually select/define the input

## *USAGE*


### *without input parameter*
        python pyfMRIqc.py
		
		
### *with input parameter*
        python pyfMRIqc.py -n <func_nift_file> -s <SNR_voxel_perc>
        python pyfMRIqc.py -n <func_nift_file> -s <SNR_voxel_perc> -m <motion_file>
        either
            python pyfMRIqc.py -n <func_nift_file> -s <SNR_voxel_perc> (-m <motion_file>) -t <mask_threshold>
        or
            python pyfMRIqc.py -n <func_nift_file> -s <SNR_voxel_perc> (-m <motion_file>) -k <mask_nift_file>
        

## *Input*
        
	-n:   functional MR nifti file
	-m:   motion parameters file of motion correction from FSL (*.par), SPM (rp*.txt) or AFNI (*.1D)
	-s:   percentage of low values outside the mask for SNR calculation
	-o:   output directory
	either
	-t:   threshold of mean values for the mask to calculate QC etc.
	or
	-k:   nifti file containing 3D binary mask to calculate QC
    
    
    
## *Output*
pyfMRIqc.py creates the following outputs:

### *Nifti files:*

	- MEAN over time
	- VAR over time
	- MASK (binary - containing voxels above the threshold or the input mask)
	- MASK4SNR (binary - lowest n percent of lowest valuues used for SNR calculation)
	- SNR signal-to-noise ratio
	- SQUARED DIFF
	- SQUARED SCALED DIFF (Squared Diff / Mean Diff)

### *a png image containing the following plots:*
	- scaled variability: Mean (over all voxel) squared diff plot over time / global mean squared diff
	- slice by slice variability: Mean (mean per slice) squared diff plot over time / global mean squared diff
	- if motion file was selected: sum of relative movements over time (z-scored)
	- scaled mean voxel intensity: mean(data/global mean) )(z-scored)
	- variance of scaled variability (z-scored)
	- min/mean/max slice variability

### *a text file:*
	-  containing an overview about scan and QC parameters

See guidelines about the output in the [Documentation](https://drmichaellindner.github.io/pyfMRIqc/)


## *Install*  
Copy the pyfMRIqc folder in a folder of your choice on your system and add the directory to your PYTHONPATH.


## *Dependencies*  

pyfMRIqc is developed in python 3.6 and the following packages need to be installed: 
numpy, nibabel, matplotloib, scipy and easygui. 
(pyfMRIqc was developed and tested with the following versions:
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
