#!/usr/bin/python

"""
pyfMRIqc.py calculates and provides information of a functional MRI nifti file for a quality check.
pyfMRIqc.py can be used by giving input parameter or without input parameter. If no input parameter are defined
pyfMRIqc.py will guide the user through input dialogs to manually select/define the input. Nevertheless, using
pyfMRIqc.py with input parameter is the recommended way of using it to gain full functionality of the tool

USAGE
    with input parameter (recommended way to use pyfMRIqc.py):
        python pyfMRIqc.py <options>

    options
        -n:   functional MR nifti file
        -s:   percentage of voxel with the lowest values (outside the mask) for SNR calculation
        either
            -t:   threshold of minimum mean values of voxel that should be included in the quality check
        or
            -k:   binary nifti mask file of voxel that should be included in the quality check
        optional
        -o:   output directory
        -m:   motion parameters file of motion correction from FSL (*.par), SPM (rp*.txt) or AFNI (*.1D).
        -x:   if -x is set the 3D and 4D nifti output files are not saved

    Example:
        python pyfMRIqc.py -n <your_functional_nii_file> -s 5 -t 200

    If no options are defined, the user is guided through input dialogs to manually specify the minimally required
    inputs

    Example:
        python pyfMRIqc.py

OUTPUT
(all output files ends with the input filename before the file extension)

png images:
    pyfMRIqc_ ... .py:
    - MEAN_<yourfile>
        showing the mean voxel intensity over time in axial slices
    - VARIANCE(_<yourfile>_thr<XXX>)
        showing variance of the voxel time courses (max threshold XXX) in axial slices
    - MASK_<yourfile>
        showing the voxels included in the QC (based on -k mask or -t threshold input) in blue and the n% of voxels with
        lowest mean voxel intensity (based on -s input) in green on top of the mean image.
    - BINMEAN_<yourfile>
        The value range of mean voxel intenstiy is devided in 50 bins with equal number of voxels.
        This image shows the average time course over voxels for each of the bins.
    - SUM_SQUARED_SCALED_DIFF_<yourfile>
    - PLOTS_<yourfile> containing the following:
        - scaled variability: Mean (over all voxel) squared diff plot over time / global mean squared diff
        - slice by slice variability: Mean (mean per slice) squared diff plot over time / global mean squared diff
        - if motion file was selected: sum of relative movements over time (z-scored)
        - scaled mean voxel intensity: mean(data/global mean) )(z-scored)
        - variance of scaled variability (z-scored)
        - min/mean/max slice variability

html file
    pyfMRIqc_HTML_<yourfile>.html containing:
    - Overview about scan and QC parameters
    - All png images explained above
    - Summary of signal to noise ratio (SNR) calculation
    - Summary of motion parameter (if motion parameter file was specified as input)

text file
    pyfMRIqc_textfile_<yourfile>.txt containing an overview about scan, QC and motion parameters (similar to the html
    file)

nifti images (if -x is not set):
    - mean_<yourfile>
        mean voxel intensity over time (3D)
    - variance_<yourfile>
        variance of the voxel time courses(3D)
    - mask_<yourfile>
        binary - containing voxels above the threshold or the input mask (3D)
    - mask4snr_<yourfile>
        binary - lowest n percent of lowest values used for SNR calculation (3D)
    - snr_<yourfile>
        voxel-wise signal-to-noise ratio (3D)
    - squared_scale_diff_<yourfile>
        squared scaled signal variability: squared difference between two consecutive volumes divided by the global
        mean difference

LICENCE
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License (GPLv3) as published
by the Free Software Foundation;

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.

AUTHORS
Michael Lindner and Brendan Williams
University of Reading, 2019
School of Psychology and Clinical Language Sciences
Center for Integrative Neuroscience and Neurodynamics
"""

import os
import sys
import getopt
import nibabel as nib
import numpy as np
import copy
from datetime import datetime
from scipy import stats
import easygui
import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt

vers = 1.1


def main():
    # check input options
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hn:m:t:k:s:o:x", ["help"])
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
    niioutput = 1

    if len(opts) > 0:

        # check parameter input
        for o, a in opts:
            if o == "-n":
                niifile = str(a)
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
            elif o == "-x":
                niioutput = 0

        if motionfile == '':
            motionfile = None
        if snrvoxelpercentage == '':
            snrvoxelpercentage = 5

        # check possible input errors and minimum requirements
        if niifile == '':
            print("INPUT ERROR (missing input): functional nifti file needs to be defined (-n). See pyfMRIqc.py -h")
            sys.exit(2)
        if maskniifile is not None and maskthresh is not None:
            print("INPUT ERROR: Only mask (-k) or threshold (-t) can be defined. See pyfMRIqc.py -h")
            sys.exit(2)
        if maskniifile is None and maskthresh is None:
            print("INPUT ERROR (missing input): mask (-k) or threshold (-t) need to be defined. See pyfMRIqc.py -h")
            sys.exit(2)

    else:
        # input dialogs if no files are give as input parameter

        # Get file extensions
        niiext = easygui.enterbox(title='Input file extension', msg='Specify file extension as either *.nii or *.nii.gz',
                                  default = "*.nii.gz")

        motionext = easygui.enterbox(title='Input file extension', msg='Optional: Specify file extension for motion file'
                                                                       ' if used as *.par for FSL, *.txt for SPM, or *.1D'
                                                                       'for AFNI', default = "*.par")

        if motionext == None:
            motionext = "*.par"

        # functional file
        if niifile == '':
            niifile = easygui.fileopenbox(title='Select functional image (*.nii or *.nii.gz)', multiple=False, default = niiext)

        # motion file
        if motionfile == '':
            motionfile = easygui.fileopenbox(title='Optional: Select motion correction files (from FSL, SPM or AFNI) - '
                                                   'Press cancel to ignore the input', multiple=False, default = motionext)

        # mask threshold
        if maskthresh is None:
            maskthresh = easygui.enterbox(title='Input mask threshold', msg='Specify threshold (mean value) for mask',
                                          default='200')

            if maskthresh is not None:
                maskthresh = int(maskthresh)

        # mask nifti file
        if maskniifile is None and maskthresh is None:
            maskniifile = easygui.fileopenbox(title='Select mask file', multiple=False, default="*.nii")

        # get SNR percentage value
        snrmsg = 'Specify as percentage of low value voxels for SNR calculation. (e.g. 25% = 25 not 0.25)'
        snrvoxelpercentage = easygui.enterbox(title='Percentage of low value voxels for SNR calculation',
                                              msg=snrmsg, default='5')
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
        outputdirectory = os.path.join(filepath, "pyfMRIqc_" + fname)

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

    process(niifile, motionfile, maskthresh, maskniifile, outputdirectory, fname, fext, snrvoxelpercentage, niioutput)


# noinspection PyBroadException
def process(niifile, motionfile, maskthresh, maskniifile, outputdirectory, fname, fext, snrvoxelpercentage, niioutput):
    # Load and get func data
    print("Load File")
    nii = nib.load(niifile)
    data = nii.get_fdata()
    shape = np.array(data)[:, :, :, 0].shape
    header = nii.header
    affine = nii.affine
    voxelsize = header['pixdim'][1]
    nrvoxel = shape[0] * shape[1] * shape[2]
    nrvolumes = np.array(data)[0, 0, 0, :].shape
    prefix = "pyfMRIqc_"
    del nii

    # Create mean data and nifti image
    print("MEAN")
    meandata = np.mean(data, axis=3)
    header2 = header
    header2['glmin'] = np.min(meandata)
    header2['glmax'] = np.max(meandata)
    new_img = nib.Nifti1Image(meandata, affine, header2)
    newfilename = os.path.join(outputdirectory, prefix + "mean_" + fname + fext)
    if niioutput == 1:
        nib.save(new_img, newfilename)

    # create png
    pngfilename = os.path.join(outputdirectory, prefix + 'MEAN_' + fname + '.png')
    meanimage = nii2image(meandata, 'Mean', pngfilename)
    meanimage = meanimage / np.max(meanimage[:]) * 255
    nrbins = 50

    # bins equal number of voxel
    binmean = np.zeros((nrbins, nrvolumes[0]))
    d = meandata.reshape((nrvoxel,))
    di = np.argsort(d)
    dbins = np.linspace(0, nrvoxel, nrbins)
    dbinlabs = np.zeros((nrbins,))

    for nn in range(nrbins):
        m = np.zeros((nrvoxel, 1))
        try:
            m[int(np.round(dbins[nn], decimals=0)):int(np.round(dbins[nn+1], decimals=0))] = 1
        except:
            m[int(np.round(dbins[nn], decimals=0)):] = 1
        binmask = m[di.argsort()]
        binmask = binmask.reshape((shape[0], shape[1], shape[2]))
        k = np.where(binmask == 1)
        if len(k[0]) > 0:
            binmean[nn, :] = stats.zscore(np.mean(data[k[0][:], k[1][:], k[2][:], :], axis=0))
            dbinlabs[nn] = np.max(meandata[k[0][:], k[1][:], k[2][:]])
        del k
        del binmask
        del m

    binticks = np.arange(1, nrbins, 5)
    binticklabels = dbinlabs[1::5]
    binticklabels = binticklabels.astype(int)

    # create image
    plt.figure(num=None, figsize=(5, 3.5), dpi=300, facecolor=(210 / 255, 227 / 255, 244 / 255), edgecolor='k')
    plt.imshow(binmean, cmap='gray', aspect='auto')
    plt.xlabel("volumes")
    plt.ylabel("bins")
    plt.yticks(binticks, binticklabels, fontsize=5)
    plt.xticks(fontsize=5)
    plt.title("bins with equal number of voxel (" + str(int(nrvoxel/nrbins)) + " per bin)")
    pngfilename = os.path.join(outputdirectory, prefix + 'BINMEAN_' + fname + '.png')
    plt.savefig(pngfilename, dpi=300, facecolor=(210 / 255, 227 / 255, 244 / 255))

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
    snrmax = vec[0][int(nrvoxel * snrvoxelpercentage / 100)]
    snrmask = np.where((meandata4snr < vec[0][int(nrvoxel * snrvoxelpercentage / 100)]) & (meandata4snr > 0), 1, 0)

    # Mean noise for SNR
    # meannoise = np.mean(meandata4snr[snrmask == 1])
    stdnoise = np.std(meandata4snr[snrmask == 1])
    # Save SNR mask
    new_img_snr = nib.Nifti1Image(snrmask, affine, header)
    newfilename = os.path.join(outputdirectory, prefix + "mask4snr_" + fname + fext)
    if niioutput == 1:
        nib.save(new_img_snr, newfilename)

    mask = []
    # Create or load mask depending on user input
    if maskthresh is not None:  # in threshold of mask input
        # check overlap between maskthresh and snr max value
        if snrmax > maskthresh:
            print("INPUT ERROR: Given mask threshold is smaller than max intensity value of given percentage of "
                  "voxel for SNR calculation.\n"
                  "SNR value range of " + str(snrvoxelpercentage) +
                  "% voxel with lowest intensity = 0 - " + str(snrmax) + "\n"
                  "Total value range = 0 - " + str(np.max(data[:])) + "\n"
                  "Specified mask threshold = " + str(maskthresh) + "\n"
                  "Either increase threshold or decrease percentage of voxel for SNR.")
            sys.exit(2)
        # Create mask
        print("Create MASK")
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

        if np.min(meandata[mask == 1]) > maskthresh:
            print("INPUT ERROR: Mask contains voxel with smaller intensity than max intensity value of given "
                  "percentage of voxel for SNR calculation.\n"
                  "SNR value range of " + str(snrvoxelpercentage) +
                  "% voxel with lowest intensity = 0 - " + str(snrmax) + "\n"
                  "Total value range = 0 - " + str(np.max(data[:])) + "\n"
                  "Minimum voxel intensity in mask = " + str(np.min(meandata[mask == 1])) + "\n"
                  "Either change mask or decrease percentage of voxel for SNR.")
            sys.exit(2)
    # Save mask
    new_img = nib.Nifti1Image(mask, affine, header)
    newfilename = os.path.join(outputdirectory, prefix + "mask_" + fname + fext)
    if niioutput == 1:
        nib.save(new_img, newfilename)
    pngfilename = os.path.join(outputdirectory, prefix + 'MASK_' + fname + '.png')
    maskimage = nii2image(mask, 'Mask', pngfilename)

    # create mask image
    imageshape = np.array(meanimage).shape
    meanimagecol = np.zeros((imageshape[0], imageshape[1], 3))
    meanimagecol[:, :, 0] = meanimage
    # add SNR mask in green
    pngfilename = os.path.join(outputdirectory, prefix + 'SNR_' + fname + '.png')
    snrimage = nii2image(snrmask, 'SNR', pngfilename)
    meanimage0 = copy.deepcopy(meanimage)
    meanimage0[snrimage > 0] = 255
    meanimagecol[:, :, 1] = meanimage0
    # add mask in blue
    meanimage2 = copy.deepcopy(meanimage)
    meanimage2[maskimage > 0] = 255
    meanimagecol[:, :, 2] = meanimage2
    pngfilename = os.path.join(outputdirectory, prefix + 'MASK_' + fname + '.png')
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
    rm = []
    relrm = []
    rmsum = []
    nrrm01 = []
    nrrm05 = []
    nrrmv = []
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

        elif motionext == ".1D":  # AFNI
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
    # calculate standard deviation
    stddata = np.std(data, axis=3)

    # signal to noise ratio
    print("SNR")
    # calculate SNR
    snrdata = np.divide(meandata, stdnoise)
    # create png
    # pngfilename = os.path.join(outputdirectory, prefix + 'SNR.png')
    # snrimage = nii2image(snrdata, 'SNR', pngfilename)
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
    newfilename = os.path.join(outputdirectory, prefix + "snr_" + fname + fext)
    if niioutput == 1:
        nib.save(new_img, newfilename)

    del meandata
    del stddata

    # variance
    print("VAR")
    vardata = np.var(data, axis=3)
    # create png
    vardata = vardata.astype(np.int32)
    varfilename = os.path.join(outputdirectory, prefix + 'VARIANCE_' + fname + '.png')
    varthresh = nii2image(vardata, 'Variance', varfilename)
    header2 = header
    header2['glmin'] = np.nanmin(vardata)
    header2['glmax'] = np.nanmax(vardata)
    new_img = nib.Nifti1Image(vardata, affine, header2)
    newfilename = os.path.join(outputdirectory, prefix + "variance_" + fname + fext)
    if niioutput == 1:
        nib.save(new_img, newfilename)
    del vardata

    # scaled squared diff
    print("SCALED SQUARED DIFF")
    diff2data = np.power(np.diff(data, axis=3), 2)
    diff2data = diff2data.astype(np.int32)
    diff2data_a = np.divide(diff2data, np.mean(diff2data))
    del diff2data
    diff2data_a = diff2data_a.astype(np.int32)
    header2 = header
    header2['glmin'] = np.min(diff2data_a)
    header2['glmax'] = np.max(diff2data_a)
    new_img = nib.Nifti1Image(diff2data_a, affine, header2)
    newfilename = os.path.join(outputdirectory, prefix + "squared_diff_scaled_" + fname + fext)
    if niioutput == 1:
        nib.save(new_img, newfilename)

    # create SCALED SQUARED DIFF image
    sumdiff2data_a = np.sum(diff2data_a, axis=3)
    pngfilename = os.path.join(outputdirectory, prefix + 'SUM_SQUARED_DIFF_SCALED_' + fname + '.png')
    mssdthresh = nii2image(sumdiff2data_a, 'DIFF', pngfilename)

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
    plt.ylabel('mean SSD')

    # plot 2
    d2 = np.mean(np.mean(diff2data_a, axis=0), axis=0)
    plt.subplot(plotnr, 1, 2)
    plt.plot(d2.T, 'x')
    plt.xlabel('Difference image number')
    plt.ylabel('slice-wise mean SSD')

    # plot 3
    d3 = np.mean(np.mean(np.mean(np.divide(data, np.mean(data)), axis=0), axis=0), axis=0)
    vardiff = np.zeros(diff2data_a.shape[3])
    for ii in range(diff2data_a.shape[3]):
        vardiff[ii] = np.var(diff2data_a[:, :, :, ii])
    plt.subplot(plotnr, 1, 3)
    plt.plot(stats.zscore(vardiff), label='variance of SSD')
    plt.plot(stats.zscore(d3), label='normalized mean voxel intensity')
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

    plt.legend(loc='lower right')  # , shadow=True, fontsize='x-large')

    # plot 4
    d4a = np.max(d2, axis=1)
    d4b = np.min(d2, axis=1)
    d4c = np.mean(d2, axis=1)
    plt.subplot(plotnr, 1, 4)
    plt.plot(d4b, label='min')
    plt.plot(d4a, label='max')
    plt.plot(d4c, label='mean')
    plt.xlabel('Slice number')
    plt.ylabel('min/mean/max slice SSD')
    plt.legend(loc='upper right')

    # close text file
    text_file.close()

    # save and show figure
    newfilename = os.path.join(outputdirectory, prefix + "PLOTS_" + fname + ".png")
    plt.savefig(newfilename, facecolor=(210 / 255, 227 / 255, 244 / 255))
    # plt.show()

    # Open html file
    print("Create html file")
    htmlfilename = os.path.join(outputdirectory, prefix + "HTML_" + fname + ".html")
    html_output = open(htmlfilename, 'w')
    # add head to html file
    html_output.write("<html><head><body style=""background-color:#d2e3f4;""><title>pyfMRIqc output</title>")
    bootstraplines = """<script src="https://stackpath.bootstrapcdn.com/bootstrap/4.1.3/js/bootstrap.min.js" 
        integrity="sha384-ChfqqxuZUCnJSK3+MXmPNIyE6ZbWh2IMqE241rYiqJxyMiZ6OW/JmZQ5stwEULTy" 
        crossorigin="anonymous"></script> <link rel="stylesheet" 
        href="https://stackpath.bootstrapcdn.com/bootstrap/4.1.3/css/bootstrap.min.css" 
        integrity="sha384-MCw98/SFnGE8fJT3GXwEOngsV7Zt27NXFoaoApmYm81iuXoPkFOJwJ8ERdknLPMO" crossorigin="anonymous">"""
    html_output.write(bootstraplines)
    
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
        h1 { padding-top: 45px; }
        h2 { padding-top: 20px; }
        h3 { padding-top: 25px; }
        body{font-size:20px;}
        </style>
        </head>
        <body>
    """
    html_output.write(style_text)

    navbar = """
            <nav class="navbar fixed-top navbar-expand-lg navbar-light bg-light">
            <a class="navbar-brand" href="#">pyfMRIqc</a>
            <div class="collapse navbar-collapse">
                <ul class="navbar-nav">
                <li class="nav-item"><a class="nav-link" href="#Input parameter">Input parameter</a></li>
                <li class="nav-item"><a class="nav-link" href="#Scan parameter">Scan parameter</a></li>
                <li class="nav-item"><a class="nav-link" href="#QC plots">QC plots</a></li>
                <li class="nav-item"><a class="nav-link" href="#BINMEAN">BINMEAN</a></li>
                <li class="nav-item"><a class="nav-link" href="#Mean">Mean</a></li>
                <li class="nav-item"><a class="nav-link" href="#Masks">Masks</a></li>
                <li class="nav-item"><a class="nav-link" href="#Variance">Variance</a></li>
                <li class="nav-item"><a class="nav-link" href="#SNR">SNR</a></li>
                <li class="nav-item"><a class="nav-link" href="#SUMDIFF">SUMDIFF</a></li>
                <li class="nav-item"><a class="nav-link" href="#Motion">Motion</a></li>
                <li class="nav-item"><a class="nav-link" href="#About">About</a></li>
                </ul>
            </div>
            </nav>"""
    html_output.write(navbar)

    # add parameter table to html file
    html_output.write("""<div id="Input parameter"> <h1>Input parameter</h1>""")
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

    html_output.write("</div><br><hr><br>")  # horizontal line

    # Add scan parameters and user input to html file
    html_output.write("""<div id="Scan parameter"> <h1>Scan parameter</h1>""")
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
                <td>""" + str(header['pixdim'][1]) + "x" + str(header['pixdim'][2]) + "x" + str(header['pixdim'][3]) +\
                        """</td>
            </tr>
            <tr>
                <th>Total Number of Voxels</th>
                <td>""" + str(data.shape[0] * data.shape[1] * data.shape[2]) + """</td>
            </tr>
        
        </table>"""
    html_output.write(acquisition_table)

    html_output.write("</div><br><hr><br>")  # horizontal line

    # add QC plot
    html_output.write("""<div id="QC plots"> <h1>QC plots</h1>""")

    plottext = """<p>The first plot shows the mean of the scaled squared difference (SSD)over all voxel in 
                    the mask. The second plot shows the mean of the (SSD) for each slice separately. The third plot 
                    shows: 1) Normalised average of the demeaned voxel intensity of each volume. 2) Normalised 
                    variance of the SSD over all voxels in the mask. 3) (only if motion file was specified as input)
                    : Normalized sum of relative movement (sum of all relative motion translations and rotations.).
                    See <a href = "https://github.com/DrMichaelLindner/pyfMRIqc/doc/index.md#qc-plots">here</a>
                    for guidelines about how to use these plots.</p>"""
    html_output.write(plottext)
    html_output.write("""<img src ="pyfMRIqc_PLOTS_""" + fname + """.png" alt="pyfMRIqc plots" class="center">""")

    html_output.write("</div><br><hr><br>")  # horizontal line

    # Add BIN x VOLUME to html file
    html_output.write("""<div id="BINMEAN"> <h1>Mean voxel time course of bins with equal number of voxels</h1>""")
    bintext = """<p>Mean voxel time course of bins.  
                        See <a href = "https://github.com/DrMichaelLindner/pyfMRIqc/doc/index.md#mean-voxel-time-course-of-bins-with-equal-number-of-voxels">here</a>
                        for guidelines about how to use this plot.</p>"""
    html_output.write(bintext)
    html_output.write("""<img src="pyfMRIqc_BINMEAN_""" + fname +
                      """.png" alt="Mean signal from functional image" class="center">""")
    html_output.write("</div><br><hr><br>")

    # Add mean data to html file
    html_output.write("""<div id="Mean"> <h1>Mean voxel intensity</h1>""")
    html_output.write("""<img src="pyfMRIqc_MEAN_""" + fname +
                      """.png" alt="Mean signal from functional image" class="center">""")
    html_output.write("<h2> Mean voxel intensity summary</h2>")
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

    html_output.write("</div><br><hr><br>")  # horizontal line

    # Add Mask to html file
    html_output.write("<div id=\"Masks\"> <h1 class=\"sub-report-title\">Masks</h1>")
    masktext = """<p>Voxels included in the masks and used for the quality check are highlighted in blue,
                voxels used for SNR calculation are highlighted in green. 
                See <a href = "https://github.com/DrMichaelLindner/pyfMRIqc/doc/index.md#masks">here</a>
                for guidelines about how to use this plot.</p>"""
    html_output.write(masktext)
    html_output.write("""<img src="pyfMRIqc_MASK_""" + fname + """.png" alt="mask image" class="center">""")
    html_output.write("<p></p>")
    html_output.write("<h2> Mask summary</h2>")
    mean_table = """
        <table style="width:50%">
            <tr>
                <th>Total Number of Voxels</th>
                <td>""" + str(nrvoxel) + """</td>
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
    html_output.write(mean_table)

    html_output.write("</div><br><hr><br>")  # horizontal line

    html_output.write("""<div id="Variance"> <h1>Variance of voxel intensity</h1>""")
    html_output.write("<p>For better visualization the image is thresholded at max " + str(varthresh) +
                      " to minimize the scaling effect of large outliers.")
    vartext = """ See <a href = "https://github.com/DrMichaelLindner/pyfMRIqc/doc/index.md#variance-of-voxel-intensity">
                here</a> for guidelines about how to use this plot.</p>"""
    html_output.write(vartext)
    head, vfilename = os.path.split(varfilename)
    html_output.write("""<img src=""" + vfilename[:-4] + """_thr""" + str(varthresh) + vfilename[-4:] +
                      """ alt="Signal variance from functional image" class="center">""")

    html_output.write("</div><br><hr><br>")  # horizontal line

    # Add SNR data to html file
    html_output.write("""<div id="SNR"> <h1>Signal to noise ratio (SNR)</h1>""")
    snrtext = """<p>See <a href = "https://github.com/DrMichaelLindner/pyfMRIqc/doc/index.md#signal-to-noise-ratio">
                       here</a> for more info.</p><p></p>"""
    html_output.write(snrtext)

    allsclicesnrtext = ""
    for ii in range(len(snrvec)):
        s = str(snrvec[ii])
        s = s.replace('[', '')
        s = s.replace(']', '')
        s = s.replace(' ', '')
        allsclicesnrtext += s + "\n"

    snr_table = """
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
    html_output.write(snr_table)

    html_output.write("</div><br><hr><br>")  # horizontal line

    # Add Mask to html file
    html_output.write(
        "<div id=\"SUMDIFF\"> <h1 class=\"sub-report-title\">Sum of squared scaled difference over time</h1>")
    html_output.write(
        "<p>For better visualization the image is Mean squared scaled difference (SSD) is thresholded at max: " + str(mssdthresh) + ".")
    ssdtext = """This image shows the sum of the difference in voxel intensity between adjacent volumes for each 
                    voxel across all adjacent volumes acquired during scanning. See <a href = 
                    "https://github.com/DrMichaelLindner/pyfMRIqc/doc/index.md#sum-of-squared-scaled-differences-over-time">
                    here</a> for more info.</p>"""
    html_output.write(ssdtext)
    html_output.write("""<img src="pyfMRIqc_SUM_SQUARED_DIFF_SCALED_""" + fname +
                      """.png" alt="sum of squared scaled difference image" class="center">""")


    html_output.write("</div><br><hr><br>")  # horizontal line

    # Add motion parameter to html file
    html_output.write("""<div id="Motion"> <h1>Motion parameter summary</h1>""")

    if motionfile is not None:
        motion_table = """
                <table style="width:70%">
                    <tr>
                        <th>Mean absolute Movement</th>
                        <td>""" + str(np.mean(rm)) + """</td>
                    </tr>
                    <tr>
                        <th>Max absolute Movement</th>
                        <td>""" + str(np.max(rm)) + """</td>
                    </tr>
                    <tr>
                        <th>Mean relative Movement</th>
                        <td>""" + str(np.mean(relrm)) + """</td>
                    </tr>
                    <tr>
                        <th>Max relative Movement</th>
                        <td>""" + str(np.max(relrm)) + """</td>
                    </tr>
                    <tr>
                        <th>Relative movements (>0.1mm)</th>
                        <td>""" + str(len(nrrm01)) + """</td>
                    </tr>
                    <tr>
                        <th><font color=#ffa500>Relative movements (>0.5mm)</font></th>
                        <td><font color=#ffa500>""" + str(len(nrrm05)) + """</font></td>
                    </tr>
                    <tr>
                        <th><font color="red">Relative movements (>voxelsize)</font></th>
                        <td><font color="red">""" + str(len(nrrmv)) + """</font></td>
                    </tr>
                </table>"""
    else:
        motion_table = """<p>No motion parameter file was provided as input.</p>"""

    html_output.write(motion_table)

    html_output.write("</div><br><hr><br>")  # horizontal line

    # get time string
    t = datetime.now()
    strg = t.strftime('%Y/%m/%d %H:%M:%S')

    html_output.write("""<div id="About"> <h1>About</h1>""")
    about_text = """
    <br>
    <ul>
		<li>Date of quality check: """ + str(strg) + """</li>
		<li>pyfMRIqc version: """ + str(vers) + """</li>
        <li>pyfMRIqc code: 
        <a href="https://github.com/DrMichaelLindner/pyfMRIqc">https://github.com/DrMichaelLindner/pyfMRIqc</a></li>
    </ul>
    <br>
    <p><font size="4"><b>Thank you for using pyfMRIqc.py!</b></font></p>
    <p><b>AUTHORS:</b><br>
        Michael Lindner and Brendan Williams<br>
        University of Reading, 2019<br>
        School of Psychology and Clinical Language Sciences<br>
        Center for Integrative Neuroscience and Neurodynamics
        </p>
    """
    html_output.write(about_text)
    html_output.write("</div><br><hr><br>")  # horizontal line

    # close html files
    html_output.write("</body></html>")
    html_output.close()

    print("DONE!")

    print("\nThank you for using this tool!")
    print("Michael Lindner and Brendan Williams")


def nii2image(img3d, cond, pngfilename):
    matdim = np.ceil(np.sqrt(img3d.shape[2]))
    slice_x = img3d.shape[0]
    slice_y = img3d.shape[0]
    img_x = int(matdim * slice_x)
    img_y = int(matdim * slice_y)
    image = np.zeros((img_x, img_y))
    cx = 0
    cy = 0
    for ii in reversed(range(img3d.shape[2])):
        f = img3d[:, :, ii]
        image[cy:cy + slice_y, cx:cx + slice_x] = np.rot90(f)

        cx += slice_x
        if cx >= img_x:
            cx = 0
            cy += slice_y

    if cond == 'Variance':
        h = np.histogram(image, bins=1000)
        thr = h[1][max(np.argwhere(h[0] > 400))]
        vmax = thr[0]
        # title = cond + ' (threshold < ' + str(vmax) + ')'
        pngfilename = pngfilename[:-4] + "_thr" + str(vmax) + pngfilename[-4:]
    elif cond == 'DIFF':
        hs, be = np.histogram(image, bins=50, density=True)
        vmax = be[5]

    else:
        vmax = img3d.max()
        # title = cond

    vmax = np.round(vmax, decimals=0)

    sizes = np.shape(image)
    height = float(sizes[0])
    width = float(sizes[1])

    fig = plt.figure()
    fig.set_size_inches(width / height, 1, forward=False)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)

    ax.imshow(image, vmax=vmax, cmap='gray')
    if cond is not 'Mask' and cond is not 'SNR':
        plt.savefig(pngfilename, dpi=height)
        plt.close()

    # plt.show()

    if cond == 'Variance':
        return vmax
    elif cond == 'DIFF':
            return vmax
    else:
        return image


def printhelp():
    helptext = """pyfMRIqc.py calculates and provides information of a functional MRI nifti file for a quality check.
    
    USAGE
        python pyfMRIqc.py <options>

    OPTIONS
        -n:   functional MR nifti file
        -s:   percentage of voxel with the lowest values (outside the mask) for SNR calculation
        either
            -t:   threshold of minimum mean values of voxel that should be included in the quality check
        or
            -k:   binary nifti mask file of voxel that should be included in the quality check
        optional
        -o:   output directory
        -m:   motion parameters file of motion correction from FSL (*.par), SPM (rp*.txt) or AFNI (*.1D).
        -x:   if -x is set the 3D and 4D nifti output files are not saved

    Example:
        python pyfMRIqc.py -n <your_functional_nii_file> -s 5 -t 200
        
    If no options are defined, the user is guided through input dialogs to manually specify the minimally required
    inputs
    
    Example:
        python pyfMRIqc.py
    """
    print(helptext)


if __name__ == "__main__":
    main()
