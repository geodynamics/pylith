#!/usr/bin/env python

"""
This script generates a plot showing a PDF of the posterior samples from catmip.
PDFs near zero can be rendered with variable transparency for visualizing purposes.
"""

import argparse
import numpy as np
import h5py
import matplotlib.pyplot as plt


def draw_plot(FILENAME_impulses, CATMIP_results, show_plot, peak_slip, output_file, with_transparency=True):
    # Impulse locations and components
    with h5py.File(FILENAME_impulses, "r") as h5:
        impulse_slip = h5['vertex_fields/slip'][:, :, :].squeeze()  # Nimpulses x Ntotalnodes x Ncomponents
        nimpulses = impulse_slip.shape[0]

    # CATMIP results output
    with open(CATMIP_results, "rb") as fin:
        coefs = np.frombuffer(fin.read(), dtype=np.float64)  # 20,000
    inversion_coefs = coefs.reshape((nimpulses, -1))  # 200 x 100

    # Histogram with transparency option
    plt.figure(figsize=(11, 6));
    for i in range(nimpulses):
        sample_values = inversion_coefs[i, :];
        if with_transparency:
            alpha = np.min([np.abs(np.median(sample_values)) / peak_slip, 1])
        else:
            alpha = 1;
        plt.hist(sample_values, alpha=alpha);
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlabel('Slip (m)', fontsize=15);
    plt.ylabel('Number of samples', fontsize=15);
    plt.title("Inversion Results: Posterior PDFs for "+str(nimpulses)+" Slip Components", fontsize=16);
    plt.savefig(output_file);
    if show_plot:
        plt.show();
    else:
        plt.close();
    return;


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    #required
    parser.add_argument("catmip_mstop", action="store", type=int, help="Catmip mstop value (which output file to plot)")
    args = parser.parse_args()
    CATMIP_results = "output/step07a_catmip-theta%d.bin"%args.catmip_mstop
    
    FILENAME_impulses = "output/step05_greensfns-fault.h5"
    
    show_plot = True  # change to False if running without graphics
    peak_slip = 1.8;  # expected peak slip for setting transparency scale
    output_file = "output/step07-pdfs.pdf"
    draw_plot(FILENAME_impulses, CATMIP_results, show_plot, peak_slip, output_file);
