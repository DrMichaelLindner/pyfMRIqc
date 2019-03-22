#!/usr/bin/python

"""
fMRI_QC.py calculates and provides information of a functional MRI nifti file for a quality check.
fMRI_QC.py can be used by giving input parameter or without input parameter. If no input parameter are defined
fMRI_QC.py will guide the user through input dialogs to manually select/define the input

USAGE
    without input parameter:
        python fmri_qc.py

    with input parameter:
        python fmri_qc.py -n <func_nift_file> -s <SNR_voxel_perc>
        python fmri_qc.py -n <func_nift_file> -s <SNR_voxel_perc> -m <motion_file>
        either
            python fmri_qc.py -n <func_nift_file> -s <SNR_voxel_perc> (-m <motion_file>) -t <mask_threshold>
        or
            python fmri_qc.py -n <func_nift_file> -s <SNR_voxel_perc> (-m <motion_file>) -k <mask_nift_file>

INPUT
    -n:   functional MR nifti file
    -m:   motion parameters file of motion correction from FSL (*.par) or SPM (rp*.txt)
    -t:   threshold of mean values for the mask to calculate SNR etc.
    -s:   percentage of low values outside the mask for SNR calculation
    -k:   mask nifti file
    -o:   output directory

OUTPUT
It creates the following nifti images as output:
    - MEAN over time
    - VAR over time
    - MASK (binary - containing voxels above the threshold or the input mask)
    - MASK4SNR (binary - lowest n percent of lowest valuues used for SNR calculation)
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
import copy

matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
from scipy import stats
import easygui


# from tkinter import messagebox, filedialog, Tk, Text, Button, mainloop


def main():
    # check input options
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hn:m:t:k:s:o:", ["help"])
    except getopt.GetoptError:
        print('INPUT ERROR:')
        printhelp()
        sys.exit(2)

    niifile = ''
    motionfile = ''
    maskniifile = None
    maskthresh = None
    outputdirectory = ''
    snrvoxelpercentage = ''

    if len(opts) > 0:

        # check parameter input
        for o, a in opts:
            if o == "-n":
                niifile = a
            elif o == "-m":
                motionfile = a
            elif o == "-t":
                maskthresh = int(a)
            elif o == "-k":
                maskniifile = a
            elif o == "-s":
                snrvoxelpercentage = int(a)
            elif o == "-o":
                outputdirectory = a
            elif o in ("-h", "--help"):
                printhelp()
                sys.exit(2)

        if motionfile == '':
            motionfile = None
        if snrvoxelpercentage == '':
            snrvoxelpercentage = 5

        # check possible input errors and minimum requirements
        if niifile == '':
            print("INPUT ERROR (missing input): functional nifti file needs to be defined (-n). See fMRI_QC.py -h")
            sys.exit(2)
        if maskniifile is not None and maskthresh is not None:
            print("INPUT ERROR: Only mask (-k) or threshold (-t) can be defined. See fMRI_QC.py -h")
            sys.exit(2)
        if maskniifile is None and maskthresh is None:
            print("INPUT ERROR (missing input): mask (-k) or threshold (-t) need to be defined. See fMRI_QC.py -h")
            sys.exit(2)

    else:
        # input dialogs if no files are give as input parameter
        # functional file
        if niifile == '':
            niifile = easygui.fileopenbox(title='Select functional image', multiple=False, default="*.nii")

        # motion file
        if motionfile == '':
            motionfile = easygui.fileopenbox(title='Select motion correction files (from FSL or SPM)', multiple=False,
                                             default="*.nii")

        # mask threshold
        if maskthresh == None:
            maskthresh = easygui.enterbox(title='Input mask threshold', msg='Specify threshold (mean value) for mask',
                                          default='200')

            if maskthresh is not None:
                maskthresh = int(maskthresh)

        # mask nifti file
        if maskniifile == None and maskthresh == None:
            maskniifile = easygui.fileopenbox(title='Select mask file', multiple=False, default="*.nii")

        # get SNR percentage value
        snrvoxelpercentage = easygui.enterbox(title='Percentage of low value voxels for SNR calculation',
                                              msg='Specify as percentage of low value voxels for SNR calculation. (e.g. 25% = 25 not 0.25)',
                                              default='5')
        snrvoxelpercentage = snrvoxelpercentage.replace('%', '')
        snrvoxelpercentage = int(snrvoxelpercentage)
        if snrvoxelpercentage is None:
            snrvoxelpercentage = 5
        elif snrvoxelpercentage == 0:
            snrvoxelpercentage = 5

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
    if maskniifile is not None:
        print('mask nifti file: ' + maskniifile)

    process(niifile, motionfile, maskthresh, maskniifile, outputdirectory, fname, fext, snrvoxelpercentage)


def process(niifile, motionfile, maskthresh, maskniifile, outputdirectory, fname, fext, snrvoxelpercentage):
    # Load and get func data
    print("Load File")
    nii = nib.load(niifile)
    data = nii.get_fdata()
    shape = np.array(data)[:, :, :, 0].shape
    header = nii.header
    affine = nii.affine
    voxelsize = header['pixdim'][1]
    nrvoxel = shape[0] * shape[1] * shape[2]
    prefix = "fMRI_QC_"
    del nii

    # Create mean data and nifti image
    print("Create and save MEAN")
    meandata = np.mean(data, axis=3)
    #meandata = meandata.astype(np.int16)
    header2 = header
    header2['glmin'] = np.min(meandata)
    header2['glmax'] = np.max(meandata)
    new_img = nib.Nifti1Image(meandata, affine, header2)
    newfilename = os.path.join(outputdirectory, prefix + "MEAN_" + fname + fext)
    nib.save(new_img, newfilename)
    # create png
    pngfilename = os.path.join(outputdirectory, prefix + 'Mean.png')
    meanimage = nii2image(meandata, 'Mean', pngfilename)

    # Create SNR mask
    meandata4snr = np.mean(data, axis=3)
    # Get lower n precent
    vec = meandata4snr.reshape(1, nrvoxel)
    vec = np.sort(vec, axis=1)
    snrvoxelmean = np.mean(vec[0][0:int(nrvoxel * snrvoxelpercentage / 100)])
    snrvoxelstd = np.std(vec[0][0:int(nrvoxel * snrvoxelpercentage / 100)])
    snrvoxelmin = vec[0][0]
    snrvoxelmax = vec[0][int(nrvoxel * snrvoxelpercentage / 100)]
    snrvoxelnr = int(nrvoxel * snrvoxelpercentage / 100)
    snrmask = np.where((meandata4snr < vec[0][int(nrvoxel * snrvoxelpercentage / 100)]) & (meandata4snr > 0), 1, 0)
    # Mean noise for SNR
    meannoise = np.mean(meandata4snr[snrmask == 1])
    # Save SNR mask
    new_img_snr = nib.Nifti1Image(snrmask, affine, header)
    newfilename = os.path.join(outputdirectory, prefix + "MASK4SNR_" + fname + fext)
    nib.save(new_img_snr, newfilename)

    # Create or load mask depending on user input
    if maskthresh is not None:  # in threshold of mask input
        # Create mask
        print("Create and save MASK")
        mask = np.where(meandata >= maskthresh, 1, 0)
    elif maskniifile is not None:  # in case of mask input
        # Load mask nifti file
        print("Load Nifti Mask")
        masknii = nib.load(maskniifile)
        maskdata = masknii.get_fdata()
        maskshape = np.array(maskdata)[:, :, :, 0].shape
        mask = np.asarray(maskdata)
        mask = mask[:, :, :, 0]
        del masknii
        # nifti mask validation
        print("Check Nifti Mask")
        maskunique = np.unique(maskdata)
        if np.min(maskunique) < 0:
            easygui.msgbox("Warning", "Mask contains value < 0")
        elif np.max(maskunique) > 1:
            easygui.msgbox("Warning", "Mask contains value > 1")
        elif len(maskunique) != 2:
            easygui.msgbox("Warning", "Mask is not binary")
        elif sum(np.asarray(shape) - np.asarray(maskshape)) != 0:
            easygui.msgbox("Mask dimensions not equal to functional data")

    # Save mask
    new_img = nib.Nifti1Image(mask, affine, header)
    newfilename = os.path.join(outputdirectory, prefix + "MASK_" + fname + fext)
    nib.save(new_img, newfilename)
    pngfilename = os.path.join(outputdirectory, prefix + 'Mask.png')
    maskimage = nii2image(mask, 'Mask', pngfilename)

    # create mask image
    imageshape = np.array(meanimage).shape
    meanimagecol = np.zeros((imageshape[0], imageshape[1], 3))

    meanimagecol[:, :, 0] = meanimage

    # add SNR mask in green
    pngfilename = os.path.join(outputdirectory, prefix + 'SNR.png')
    snrimage = nii2image(snrmask, 'SNR', pngfilename)
    meanimage0 = copy.deepcopy(meanimage)
    meanimage0[snrimage > 0] = 255
    meanimagecol[:, :, 1] = meanimage0

    # add mask in blue
    meanimage2 = copy.deepcopy(meanimage)
    meanimage2[maskimage > 0] = 255
    meanimagecol[:, :, 2] = meanimage2
    # plt.imshow(meanimagecol / 255.0)
    pngfilename = os.path.join(outputdirectory, prefix + 'Mask.png')
    # plt.savefig(pngfilename)

    sizes = np.shape(meanimagecol)
    height = float(sizes[0])
    width = float(sizes[1])
    fig = plt.figure()
    fig.set_size_inches(width / height, 1, forward=False)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax.imshow(meanimagecol / 255.0)
    plt.savefig(pngfilename, dpi=height)


    # Apply mask
    s = np.asarray(np.asarray(np.array(data)).shape)
    # data = np.empty((s[0], s[1], s[2], s[3]))
    for ii in range(s[3]):
        data[:, :, :, ii] = np.multiply(mask, data[:, :, :, ii])

    # Open text file
    print("Create text file")
    textfilename = os.path.join(outputdirectory, prefix + "textfile_" + fname + ".txt")
    text_file = open(textfilename, "w")
    # Add scanparameter and user input text file
    text_file.write("SCAN PARAMETERS AND USER INPUT: \n\n")
    text_file.write("Slice Resolution: " + str(header['dim'][1]) + "x" + str(header['dim'][2]) + "\n")
    text_file.write("Num slices: " + str(header['dim'][3]) + "\n")
    text_file.write("Num Volumes: " + str(data.shape[3]) + "\n")
    text_file.write(
        "Voxel Size: " + str(header['pixdim'][1]) + "x" + str(header['pixdim'][2]) + "x" + str(header['pixdim'][3])
        + "\n")
    text_file.write("Total Num Voxels: " + str(data.shape[0] * data.shape[1] * data.shape[2]) + "\n")
    text_file.write("Mask Threshold: " + str(maskthresh) + "\n")
    text_file.write("Num Mask Voxels: " + str(np.sum(mask)) + "\n")
    text_file.write("Percentage of voxel with lowest values for SNR: " + str(snrvoxelpercentage) + "\n")
    text_file.write("  SNR nr of lowest voxel: " + str(snrvoxelnr) + "\n")
    text_file.write("  SNR voxel MEAN: " + str(snrvoxelmean) + "\n")
    text_file.write("  SNR voxel STD: " + str(snrvoxelstd) + "\n")
    text_file.write("  SNR voxel value range: " + str(snrvoxelmin) + " - " + str(snrvoxelmax) + "\n")

    if maskniifile is not None:
        text_file.write("Mask file: " + maskniifile + "\n")

    text_file.write("\n------------------------------------ \n\n")
    text_file.write("QC OVERVIEW: \n\n")
    text_file.write("Mean: " + str(np.mean(data)) + "\n")
    text_file.write("Mean (mask): " + str(np.mean(data[mask == 1])) + "\n")
    text_file.write("SD: " + str(np.std(data)) + "\n")
    text_file.write("SD (mask): " + str(np.std(data[mask == 1])) + "\n")

    # Load and get motion data if specified
    if motionfile is not None:
        print("Load Motion Parameters")
        # Load motion file
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
    # calculate SNR
    snrdata = np.divide(meandata, meannoise)
    # create png
    #pngfilename = os.path.join(outputdirectory, prefix + 'SNR.png')
    #snrimage = nii2image(snrdata, 'SNR', pngfilename)
    # mean SNR over slice
    snrvec = np.zeros((snrdata.shape[2], 1))
    for ii in range(snrdata.shape[2]):
        snrslice = snrdata[:, :, ii]
        maskslice = mask[:, :, ii]
        if np.sum(maskslice) > 0:
            snrvec[ii] = np.nanmean(snrslice[maskslice == 1])
        else:
            snrvec[ii] = np.nan

    # add snr to text file
    text_file.write("\nSignal-To-Noise Ratio: \n")
    text_file.write("Min Slice SNR: " + str(np.nanmin(snrvec)) + "\n")
    text_file.write("Max Slice SNR: " + str(np.nanmax(snrvec)) + "\n")
    text_file.write("ALL Slice SNRs: ")
    for ii in range(len(snrvec)):
        text_file.write(str(snrvec[ii]) + " ")
    text_file.write("\n")
    text_file.write("Mean voxel SNR: " + str(np.nanmean(snrvec)) + "\n")

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
    # create png
    vardata = vardata.astype(np.int32)
    varfilename = os.path.join(outputdirectory, prefix + 'Variance.png')
    varthresh = nii2image(vardata, 'Variance', varfilename)
    header2 = header
    header2['glmin'] = np.nanmin(vardata)
    header2['glmax'] = np.nanmax(vardata)
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
    plotnr = 4
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

    # Open html file
    print("Create html file")
    htmlfilename = os.path.join(outputdirectory, prefix + "html_" + fname + ".html")
    html_output = open(htmlfilename, 'w')
    # add head to html file
    html_output.write("<html><head><body style=""background-color:#ddd;""><title>fMRI_QC output</title></head>")

    # add title to html file
    html_output.write("<body> <h1> fMRI_QC output </h1> ")

    # add style for tables to html file
    style_text = """
        <style>
        table, th, td {
            border: 1px solid black;
            border-collapse: collapse;
            width: 400;    
            background-color: #fff;
            align="center"
        }
        th, td {
            padding: 15px;
            text-align: left;
        }
        table th {
            background-color: #eee;
        }
        hr {
            display: block;
            height: 1px;
            border: 20;
            border-top: 1px 
            solid #333;
            margin: 1em 0;
            padding: 0; 
        }
        .center {
            display: block;
            margin-left: auto;
            margin-right: auto;
            width: 90%;
}
        </style>
    """
    html_output.write(style_text)

    # add parameter table to html file
    html_output.write("<h2> fMRI_QC Input parameters </h2>")
    parameter_table = """
        <table>
            <tr>
                <th>Functional file</th>
                <td>""" + fname + """</td>
            </tr>
            <tr>
                <th>Motion file</th>
                <td>""" + str(motionfile) + """</td>
            </tr>
            <tr>
                <th>Threshold value</th>
                <td>""" + str(maskthresh) + """</td>
            </tr>
            <tr>
                <th>Mask file</th>
                <td>""" + str(maskniifile) + """</td>
            </tr>
            <tr>
                <th>SNR threshold</th>
                <td>""" + str(snrvoxelpercentage) + """</td>
            </tr>
        </table>"""
    html_output.write(parameter_table)

    html_output.write("\n<hr>\n") # horizontal line

    # Add scan parameters and user input to html file
    html_output.write("<h2> Functional image parameters </h2>")
    acquisition_table = """
        <table>
            <tr>
                <th>Slice Dimensions</th>
                <td>""" + str(header['dim'][1]) + "x" + str(header['dim'][2]) + """</td>
            </tr>
            <tr>
                <th>Number of Slices</th>
                <td>""" + str(header['dim'][3]) + """</td>
            </tr>
            <tr>
                <th>Number of Volumes</th>
                <td>""" + str(data.shape[3]) + """</td>
            </tr>
            <tr>
                <th>Voxel Size</th>
                <td>""" + str(header['pixdim'][1]) + "x" + str(header['pixdim'][2]) + "x" + str(header['pixdim'][3]) + """</td>
            </tr>
            <tr>
                <th>Total Number of Voxels</th>
                <td>""" + str(data.shape[0] * data.shape[1] * data.shape[2]) + """</td>
            </tr>
        
        </table>"""
    html_output.write(acquisition_table)

    html_output.write("\n<hr>\n")  # horizontal line

    # add QC plot
    html_output.write("<h2> Quality check plots </h2>")
    html_output.write("""<img src ="fMRI_QC_plots_""" + fname + """.png" alt="fMRI_QC plots" class="center">""")

    html_output.write("\n<hr>\n")  # horizontal line

    # Add mean data to html file
    html_output.write("<h2> Mean voxel intensity </h2>")
    html_output.write("""<img src="fMRI_QC_Mean.png" alt="Mean signal from functional image" class="center">""")
    html_output.write("<h3> Mean voxel intensity summary</h3>")
    mean_signal_table = """
        <table style="width:50%">
            <tr>
                <th>Mean Signal (unmasked)</th>
                <td>""" + str(np.mean(data)) + """</td>
            </tr>
            <tr>
                <th>Mean Signal SD (unmasked)</th>
                <td>""" + str(np.std(data)) + """</td>
            </tr>
            <tr>
                <th>Mean Signal (masked)</th>
                <td>""" + str(np.mean(data[mask == 1])) + """</td>
            </tr>
            <tr>
                <th>Mean Signal SD (masked)</th>
                <td>""" + str(np.std(data[mask == 1])) + """</td>
            </tr>
        </table>"""
    html_output.write(mean_signal_table)

    html_output.write("\n<hr>\n")  # horizontal line

    # Add Mask to html file
    html_output.write("<h2> Masks </h2>")
    html_output.write(
        "<p>Voxels in cluded inthe masks are shown in blue, voxels used for SNR calcualtion are shown in green:</p>")
    html_output.write("""<img src="fMRI_QC_MASK.png" alt="mask image" class="center">""")
    html_output.write("<p></p>")

    MEAN_table = """
        <table style="width:50%">
            <tr>
                <th>Total Number of Voxels</th>
                <td>""" + str(np.sum(mask)) + """</td>
            </tr>
            <tr>
                <th>Mask Threshold Value</th>
                <td>""" + str(maskthresh) + """</td>
            </tr>
            <tr>
                <th>Number of Masked Voxels</th>
                <td>""" + str(np.sum(mask)) + """</td>
            </tr>
        </table>"""
    html_output.write(MEAN_table)


    html_output.write("\n<hr>\n")  # horizontal line

    html_output.write("<h2> Voxel variance </h2>")
    html_output.write("Voxel variance is thresholded at max " + str(varthresh) + ":<p></p>")
    head, vfilename = os.path.split(varfilename)
    html_output.write("""<img src=""" + vfilename[:-4] + """_thr""" + str(varthresh) + vfilename[-4:] +
                      """ alt="Signal variance from functional image" class="center">""")

    html_output.write("\n<hr>\n")  # horizontal line

    # Add SNR data to html file
    html_output.write("<h2> Signal to noise ratio </h2>")
    # html_output.write("""<img src="fMRI_QC_SNR.png" alt="SNR from functional image" class="center">""")
    # html_output.write("<h3> Voxel signal to noise ratio summary</h3>")
    allsclicesnrtext = ""
    for ii in range(len(snrvec)):
        s = str(snrvec[ii])
        s = s.replace('[', '')
        s = s.replace(']', '')
        s = s.replace(' ', '')
        allsclicesnrtext += s + "\n"

    SNR_table = """
        <table style="width:70%">
            <tr>
                <th>Mask Threshold Value</th>
                <td>""" + str(maskthresh) + """</td>
            </tr>
            <tr>
                <th>Number of Masked Voxels</th>
                <td>""" + str(np.sum(mask)) + """</td>
            </tr>
            <tr>
                <th>SNR Threshold</th>
                <td>""" + str(snrvoxelpercentage) + """</td>
            </tr>
            <tr>
                <th>Number of voxels below threshold for SNR</th>
                <td>""" + str(snrvoxelnr) + """</td>
            </tr>
            <tr>
                <th>Mean values of voxel for SNR</th>
                <td>""" + str(snrvoxelmean) + """</td>
            </tr>
            <tr>
                <th>STD of voxel for SNR</th>
                <td>""" + str(snrvoxelstd) + """</td>
            </tr>
            <tr>
                <th>Value range of voxel for SNR</th>
                <td>""" + str(snrvoxelmin) + " - " + str(snrvoxelmax) + """</td>
            </tr>
            <tr>
                <th>Mean voxel SNR</th>
                <td>""" + str(np.nanmean(snrvec)) + """</td>
            </tr>
            <tr>
                <th>Min Slice SNR</th>
                <td>""" + str(np.nanmin(snrvec)) + """</td>
            </tr>
            <tr>
                <th>Min Slice SNR</th>
                <td>""" + str(np.nanmax(snrvec)) + """</td>
            </tr>
            <tr>
                <th>ALL Slice SNR</th>
                <td>""" + str(allsclicesnrtext) + """</td>
            </tr>
        </table>"""
    html_output.write(SNR_table)





    # close html files
    html_output.write("</body></html>")
    html_output.close()

    print("DONE!")

    print("\nThank you for using this tool!")
    print("    Michael Lindner")


def nii2image(img3D, cond, pngfilename):
    matdim = np.ceil(np.sqrt(img3D.shape[2]))
    slice_x = img3D.shape[0]
    slice_y = img3D.shape[0]
    img_x = int(matdim * slice_x)
    img_y = int(matdim * slice_y)
    image = np.zeros((img_x, img_y))
    cx = 0
    cy = 0
    for ii in reversed(range(img3D.shape[2])):
        f = img3D[:, :, ii]
        image[cy:cy + slice_y, cx:cx + slice_x] = np.rot90(f)

        cx += slice_x
        if cx >= img_x:
            cx = 0
            cy += slice_y

    if cond == 'Variance':
        h = np.histogram(image, bins=1000)
        thr = h[1][max(np.argwhere(h[0] > 400))]
        vmax = thr[0]
        title = cond + ' (threshold < ' + str(vmax) + ')'
        pngfilename = pngfilename[:-4] + "_thr" + str(vmax) + pngfilename[-4:]
    else:
        vmax = img3D.max()
        title = cond

    sizes = np.shape(image)
    height = float(sizes[0])
    width = float(sizes[1])

    fig = plt.figure()
    fig.set_size_inches(width / height, 1, forward=False)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)

    #dpi = 300
    #margin = 0.05
    #figsize = (1 + margin) * img_y / dpi, (1 + margin) * img_x / dpi
    #fig = plt.figure(figsize=figsize, dpi=dpi)
    #extent = (0, img_x * 3, img_y * 3, 0)

    #plt.imshow(image, vmax=vmax, extent=extent, interpolation=None, cmap='gray')
    #plt.axis('off')
    #plt.title(title)
    ax.imshow(image, vmax=vmax, cmap='gray')
    if cond is not 'Mask' and cond is not 'SNR':
        plt.savefig(pngfilename, dpi=height)
        plt.close()

    # plt.show()

    if cond == 'Variance':
        return (vmax)
    else:
        return (image)


def printhelp():
    helptext = """fMRI_QC.py calculates and provides information of a functional MRI nifti file for a quality check.
        fMRI_QC.py can be used by giving input parameter or without input parameter. If no input parameter are defined 
        fMRI_QC.py will guide the user through input dialogs to manually select/define the input 

    USAGE
        without input parameter:
            python fmri_qc.py

        with input parameter:
            python fmri_qc.py -n <func_nift_file> -s <SNR_voxel_perc>
            python fmri_qc.py -n <func_nift_file> -s <SNR_voxel_perc> -m <motion_file>
            either
                python fmri_qc.py -n <func_nift_file> -s <SNR_voxel_perc> (-m <motion_file>) -t <mask_threshold>
            or
                python fmri_qc.py -n <func_nift_file> -s <SNR_voxel_perc> (-m <motion_file>) -k <mask_nift_file>

    INPUT
        -n:   functional MR nifti file
        -m:   motion parameters file of motion correction from FSL (*.par) or SPM (rp*.txt)
        -s:   percentage of low values outside the mask for SNR calculation
        -t:   threshold of mean values for the mask to calculate SNR etc.
        -k:   mask nifti file
        -o:   output directory


    OUTPUT
    It creates the following nifti images as output:
        - MEAN over time
        - VAR over time
        - MASK (binary - containing voxels above the threshold or the input mask)
        - MASK4SNR (binary - lowest n percent of lowest valuues used for SNR calculation)
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
