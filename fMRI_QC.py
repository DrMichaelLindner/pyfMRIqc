#!/usr/bin/python

"""
fMRI_QC.py calculates and provides information of a functional MRI nifti file for a quality check.

USAGE
    python fmri_qc.py
    python fmri_qc.py -n <func_nift_file>
    python fmri_qc.py -n <func_nift_file> -m <motion_file>
    python fmri_qc.py -n <func_nift_file> -m <motion_file> -t <mask_threshold>
    python fmri_qc.py -n <func_nift_file> -m <motion_file> -t <mask_threshold> -k <mask_nift_file>
    python fmri_qc.py -n <func_nift_file> -m <motion_file> -t <mask_threshold> -k <mask_nift_file> -o <output_path>

INPUT
    -n:   functional MR nifti file
    -m:   motion parameters file of motion correction from FSL (*.par) or SPM (rp*.txt)
    -t:   threshold of mean values for the mask to calculate SNR etc.
    -k:   mask nifti file
    -o:   output folder
All input are optionally: if not defined input dialogs will pop up to select the files manually.

OUTPUT
It creates the following nifti images as output:
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
    - containing an overview about scan and QC parameters

LICENCE
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License (GPLv3) as published
by the Free Software Foundation;

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.

AUTHOR
Michael Lindner
University of Reading, 2018
School of Psychology and Clinical Language Sciences
Center for Integrative Neuroscience and Neurodynamics
"""

import os
import sys
import getopt
import nibabel as nib
import numpy as np
import matplotlib

matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
from scipy import stats
import easygui
# from tkinter import messagebox, filedialog, Tk, Text, Button, mainloop


def main():
    # check input options
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hn:m:t:k:o:", ["help"])
    except getopt.GetoptError:
        print('python fmri_qc.py -n <func_nift_file> -m <motion_file> -t <mask_threshold> -k <mask_nift_file>')
        sys.exit(2)

    niifile = ''
    motionfile = ''
    masknii = ''
    maskthresh = ''
    outputdirectory = ''

    for o, a in opts:
        if o == "-n":
            niifile = a
        elif o == "-m":
            motionfile = a
        elif o == "-t":
            maskthresh = int(a)
        elif o == "-k":
            masknii == a
        elif o == "-o":
            outputdirectory = a
        elif o in ("-h", "--help"):
            printhelp()
            sys.exit(2)

    # input dialogs if no files are give as input parameter
    # functional file
    if niifile == '':
        niifile = easygui.fileopenbox(title='Select functional image', multiple=False, default="*.nii")
        # niifile = filedialog.askopenfilename(title="Select functional image", multiple=False,
        #                              filetypes=(("nifti", "*.nii"),("all files","*.*")))
    # motion file
    if motionfile == '':
        motionfile = easygui.fileopenbox(title='Select motion correction files (from FSL or SPM)', multiple=False,
                                         default="*.nii")
        # motionfile = filedialog.askopenfilename(title="Select motion correction files (from FSL or SPM)", multiple=False,
        #                                      filetypes=(("nifti", "*.nii"), ("all files", "*.*")))

    # mask threshold
    if maskthresh == '':
        maskthresh = easygui.integerbox(title='Input mask threshold', msg='Specify threshold (mean value) for mask',
                                      default='200')

        # maskthresh = int(maskthresh)
        # thresh = Tk()
        #
        # def retrieve_input():
        #     inputValue = textBox.get("1.0", "end-1c")
        #     maskthresh = int(inputValue)
        #     print(maskthresh)
        #
        # textBox = Text(thresh, height=1, width=25)
        # textBox.pack()
        # buttonCommit = Button(thresh, height=1, width=20, text="Input mask threshold",
        #                       command=lambda: retrieve_input())
        # buttonCommit.pack()
        # Button(thresh, height=1, width=20, text="Close input window", command=thresh.destroy).pack()
        # mainloop()

    # mask nifti file
    if masknii == '' and maskthresh == None:
       masknii = easygui.fileopenbox(title='Select mask file', multiple=False, default="*.nii")
       maskthresh = 0
       # motionfile = filedialog.askopenfilename(title="Select mask file", multiple=False,
       #                                      filetypes=(("nifti", "*.nii"), ("all files", "*.*")))

    # get filename
    filepath, filename = os.path.split(niifile)
    x, fx = os.path.splitext(filename)
    if fx == ".gz":
        fname, fx2 = os.path.splitext(x)
        fext = ".nii.gz"
    else:
        fname = x
        fext = ".nii"

    # output folder
    if outputdirectory == '':
        filepath, filename = os.path.split(niifile)
        outputdirectory = os.path.join(filepath, "fMRI_QC_" + fname)

    # create output folder if not existing
    if not os.path.exists(outputdirectory):
        os.makedirs(outputdirectory)

    # print input
    print('Functional MRI nifti file: ' + niifile)
    if motionfile is not None:
        print('Motion parameter file: ' + motionfile)
    print('Mask threshold: ' + str(maskthresh))
    if masknii is not None:
        print('mask nifti file: ' + masknii)

    process(niifile, motionfile, maskthresh, masknii, outputdirectory, fname, fext)


def process(niifile, motionfile, maskthresh, masknii, outputdirectory, fname, fext):

    # load and get func data
    print("Load File")
    nii = nib.load(niifile)
    data = nii.get_fdata()
    shape = np.array(data)[:, :, :, 0].shape
    header = nii.header
    affine = nii.affine
    voxelsize = header['pixdim'][1]
    prefix = "fMRI_QC_"
    del nii

    if masknii is not None:

        # load mask nifti file
        print("Load Nifti Mask")
        mknii = nib.load(masknii)
        mkdata = mknii.get_fdata()
        mkshape = np.array(mkdata)[:, :, :, 0].shape
        del mknii

        # nifti mask validation
        mkunique = np.unique(mkdata)
        if np.min(mkunique) < 0:
            easygui.msgbox("Warning", "Mask contains value < 0")
        elif np.max(mkunique) > 1:
            easygui.msgbox("Warning", "Mask contains value > 1")
        elif len(mkunique) != 2:
            easygui.msgbox("Warning", "Mask is not binary")
        elif sum(np.asarray(shape) - np.asarray(mkshape)) != 0:
            easygui.msgbox("Mask dimensions not equal to functional data")

        # Apply mask

        s = np.asarray(np.asarray(np.array(data)).shape)
        #data = np.empty((s[0], s[1], s[2], s[3]))
        for ii in range(s[3]):
            data[:, :, :, ii] = np.multiply(mkdata[:, :, :, 0], data[:, :, :, ii])


    # create mean data and nifti image
    print("Create and save MEAN")
    meandata = np.mean(data, axis=3)
    meandata = meandata.astype(np.int16)
    header2 = header
    header2['glmin'] = np.min(meandata)
    header2['glmax'] = np.max(meandata)
    new_img = nib.Nifti1Image(meandata, affine, header2)
    newfilename = os.path.join(outputdirectory, prefix + "MEAN_" + fname + fext)
    nib.save(new_img, newfilename)

    # create and save mask
    print("Create and save MASK")
    mask = np.where(meandata >= maskthresh, 1, 0)
    mask2 = np.where((meandata < maskthresh / 5) & (meandata > 0), 1, 0)
    new_img = nib.Nifti1Image(mask, affine, header)
    newfilename = os.path.join(outputdirectory, prefix + "MASK_" + fname + fext)
    nib.save(new_img, newfilename)

    # open text file
    print("Create text file")
    textfilename = os.path.join(outputdirectory, prefix + "textfile_" + fname + ".txt")
    text_file = open(textfilename, "w")
    text_file.write("SCAN PARAMETERS: \n\n")
    text_file.write("Slice Resolution: " + str(header['dim'][1]) + "x" + str(header['dim'][2]) + "\n")
    text_file.write("Num slices: " + str(header['dim'][3]) + "\n")
    text_file.write("Num Volumes: " + str(data.shape[3]) + "\n")
    text_file.write(
        "Voxel Size: " + str(header['pixdim'][1]) + "x" + str(header['pixdim'][2]) + "x" + str(header['pixdim'][3])
        + "\n")
    text_file.write("Total Num Voxels: " + str(data.shape[0] * data.shape[1] * data.shape[2]) + "\n")
    text_file.write("Mask Threshold: " + str(maskthresh) + "\n")
    text_file.write("Num Mask Voxels: " + str(np.sum(mask)) + "\n")

    text_file.write("------------------------------------ \n\n")
    text_file.write("QC OVERVIEW: \n\n")

    text_file.write("Mean: " + str(np.mean(data)) + "\n")
    text_file.write("Mean (mask): " + str(np.mean(data[mask == 1])) + "\n")
    text_file.write("SD: " + str(np.std(data)) + "\n")
    text_file.write("SD (mask): " + str(np.std(data[mask == 1])) + "\n")

    # load and get motion data if specified
    plotnr = 4
    if motionfile is not None:
        print("Load Motion Parameters")
        # load motion file
        x, motionext = os.path.splitext(motionfile)
        rm = np.loadtxt(motionfile)

        # degree in mm assuming head radius is 5cm
        if motionext == ".txt":  # SPM
            for ii in [3, 4, 5]:
                rm[:, ii] = rm[:, ii] * 50

        elif motionext == ".par":  # FSL
            for ii in [0, 1, 2]:
                rm[:, ii] = rm[:, ii] * 50

        # absolute values
        rm = np.absolute(rm)

        # sum absolute values
        rmsum = np.sum(rm, axis=1)

        # create relative values
        relrm = np.diff(rm, axis=1)

        # get thrtesholds:
        nrrm01 = rm[np.where(rm > .1)]
        nrrm05 = rm[np.where(rm > .5)]
        nrrmv = rm[np.where(rm > voxelsize)]

    # create nifti files
    # --------------------
    print("Create NIFTI files")
    # calculate standard deviation
    stddata = np.std(data, axis=3)

    # signal to noise ratio
    print("- SNR")
    meannoise = np.mean(meandata[mask2 == 1])
    snrdata = np.divide(meandata, meannoise)

    snrvec = np.zeros((snrdata.shape[2], 1))
    for ii in range(snrdata.shape[2]):
        snrslice = snrdata[:, :, ii]
        snrvec[ii] = np.mean(snrslice[mask[:, :, ii] == 1])

    # add snr to text file
    text_file.write("\nSignal-To-Noise Ratio: \n")
    text_file.write("Min Slice SNR: " + str(np.min(snrvec)) + "\n")
    text_file.write("Max Slice SNR: " + str(np.max(snrvec)) + "\n")
    text_file.write("ALL Slice SNRs: ")
    for ii in range(len(snrvec)):
        text_file.write(str(snrvec[ii]) + " ")
    text_file.write("\n")
    text_file.write("Mean voxel SNR: " + str(np.mean(snrvec)) + "\n")

    # save SNR nifti
    snrdata = snrdata.astype(np.int16)
    header2 = header
    header2['glmin'] = np.min(snrdata)
    header2['glmax'] = np.max(snrdata)

    new_img = nib.Nifti1Image(snrdata, affine, header2)
    newfilename = os.path.join(outputdirectory, prefix + "SNR_" + fname + fext)
    nib.save(new_img, newfilename)

    del meandata
    del stddata

    # variance
    print("- VAR")
    vardata = np.var(data, axis=3)
    vardata = vardata.astype(np.int32)
    header2 = header
    header2['glmin'] = np.min(vardata)
    header2['glmax'] = np.max(vardata)
    new_img = nib.Nifti1Image(vardata, affine, header2)
    newfilename = os.path.join(outputdirectory, prefix + "VAR_" + fname + fext)
    nib.save(new_img, newfilename)
    del vardata

    # squared diff
    print("- SQUARED DIFF")
    diff2data = np.power(np.diff(data, axis=3), 2)
    diff2data = diff2data.astype(np.int32)
    header2 = header
    header2['glmin'] = np.min(diff2data)
    header2['glmax'] = np.max(diff2data)
    new_img = nib.Nifti1Image(diff2data, affine, header2)
    newfilename = os.path.join(outputdirectory, prefix + "squared_diff_" + fname + fext)
    nib.save(new_img, newfilename)
    # del diffdata

    print("- SCALED SQUARED DIFF")
    diff2data_a = np.divide(diff2data, np.mean(diff2data))
    del diff2data
    diff2data_a = diff2data_a.astype(np.int32)
    header2 = header
    header2['glmin'] = np.min(diff2data_a)
    header2['glmax'] = np.max(diff2data_a)
    new_img = nib.Nifti1Image(diff2data_a, affine, header2)
    newfilename = os.path.join(outputdirectory, prefix + "squared_diff_scaled_" + fname + fext)
    nib.save(new_img, newfilename)

    # plot data
    # ------------------
    print("Create Plot")
    fig = plt.figure(frameon=False)
    fig.set_size_inches(12, 14)

    # plot 1
    d1 = np.mean(np.mean(np.mean(diff2data_a, axis=0), axis=0), axis=0)
    plt.subplot(plotnr, 1, 1)
    plt.plot(d1)
    plt.xlabel('Difference image number')
    plt.ylabel('scaled variability')

    # plot 2
    d2 = np.mean(np.mean(diff2data_a, axis=0), axis=0)
    plt.subplot(plotnr, 1, 2)
    plt.plot(d2.T, 'x')
    plt.xlabel('Difference image number')
    plt.ylabel('slice by slice variability')

    # plot 3
    d3 = np.mean(np.mean(np.mean(np.divide(data, np.mean(data)), axis=0), axis=0), axis=0)
    vardiff = np.zeros(diff2data_a.shape[3])
    for ii in range(diff2data_a.shape[3]):
        vardiff[ii] = np.var(diff2data_a[:, :, :, ii])
    plt.subplot(plotnr, 1, 3)
    plt.plot(stats.zscore(vardiff), label='variance of scaled variability')
    plt.plot(stats.zscore(d3), label='scaled mean voxel intensity')
    plt.xlabel('Image number')
    plt.ylabel('Normalized Amplitudes')
    if motionfile is not None:
        # add movent to textfile
        text_file.write("\nMovement: \n")
        text_file.write("Mean absolute Movement: " + str(np.mean(rm)) + "\n")
        text_file.write("Max absolute Movement: " + str(np.max(rm)) + "\n")
        text_file.write("Mean relative Movement: " + str(np.mean(relrm)) + "\n")
        text_file.write("Max relative Movement: " + str(np.max(relrm)) + "\n")
        text_file.write("Movements (>0.1mm): " + str(len(nrrm01)) + "\n")
        text_file.write("Movements (>0.5mm): " + str(len(nrrm05)) + "\n")
        text_file.write("Movements (>voxelsize): " + str(len(nrrmv)) + "\n")

        # add line to plot
        plt.plot(stats.zscore(rmsum), label='sum of relative movements')

    plt.legend(loc='upper right')  # , shadow=True, fontsize='x-large')

    # plot 4
    d4a = np.max(d2, axis=1)
    d4b = np.min(d2, axis=1)
    d4c = np.mean(d2, axis=1)
    plt.subplot(plotnr, 1, 4)
    plt.plot(d4b, label='min')
    plt.plot(d4a, label='max')
    plt.plot(d4c, label='mean')
    plt.xlabel('Slice number')
    plt.ylabel('min/mean/max slice variability')
    plt.legend(loc='upper right')

    # close text file
    text_file.close()

    # save and show figure
    newfilename = os.path.join(outputdirectory, prefix + "plots_" + fname + ".png")
    plt.savefig(newfilename)
    # plt.show()

    print("DONE!")

    print("\nThank you for using this tool!")
    print("    Michael Lindner")


def printhelp():
    helptext = """fMRI_QC.py calculates and provides information of a functional MRI nifti file for a quality check.

    USAGE
        python fmri_qc.py 
        python fmri_qc.py -n <func_nift_file> 
        python fmri_qc.py -n <func_nift_file> -m <motion_file> 
        python fmri_qc.py -n <func_nift_file> -m <motion_file> -t <mask_threshold>
        python fmri_qc.py -n <func_nift_file> -m <motion_file> -t <mask_threshold> -o <output_folder>

    INPUT
        -n:   functional MR nifti file 
        -m:   motion parameters file of motion correction from FSL (*.par) or SPM (rp*.txt)
        -t:   threshold of mean values for the mask to calculate SNR etc.
        -o:   output folder  
    All input are optionally: if not defined input dialogs will pop up to select the files manually.

    OUTPUT
    It creates the following nifti images as output:
        - MEAN over time
        - VAR over time
        - MASK (containing voxels above the threshold)
        - SNR: signal-to-noise ration
        - SQUARED DIFF (squared difference between two consecutive volumes)
        - SQUARED SCALED DIFF (Squared Diff / Mean Diff)  

    and the following plots (and saved png):
        - scaled variabiliy: Mean (over all voxel) squared diff plot over time / global mean squared diff
        - slice by slice variabiliy: Mean (per slice) squared diff plot over time / global mean squared diff
        - if motion file was selected: mean absolute motion over time
        - scaled mean voxel intensity: mean(data/global mean)
        - variance of scaled variability
        - min/mean/max slice variability
    """
    print(helptext)


if __name__ == "__main__":
    main()
