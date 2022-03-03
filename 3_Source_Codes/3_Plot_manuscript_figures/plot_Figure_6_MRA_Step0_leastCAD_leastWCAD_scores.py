#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
from datetime import datetime, timezone
from csv import reader
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def convert_str_to_datetime(date, dateformat="%d/%m/%Y", tzinfo=timezone.utc):
    refdate = datetime.strptime(date, dateformat)
    refdate = refdate.replace(tzinfo=tzinfo)
    return refdate


def read_CAD_WCAD(source='source.csv'):
    '''Method to extract data from either the CAD.csv and WCAD.csv result file
    of the Predict_SARSCOV2_COVID19_VIA_CPTPP class and return these
    data as a dictionary.

    Dictionary keys and values (note: each value  in values is a list obj):
      'batch':[xxx] - give the batch number (integer) that gave the lowest
                      CAD or WCAD score.
      'least_score':[iii, iii] - gives the iteration number of that batch and
                                 the CAD/WCAD score from that iteration
      'final_mus':[ii,ii,...] - daily mean CCPs of the iteration with the
                                lowest CAD/WCAD score.
      'inf_days_bins':[i,i,...] - the day numbers of the SARS-CoV-2 infection
                                  trend curve with the least CAD/WCAD score.
      'inf_days_ncases':[i,i,...] - the daily number of local SARS-CoV-2
                                    infectees of the iteration with the lowest
                                    CAD/WCAD score.
      'covid_days_ncases':[i,i,...] - the daily number of local COVID-19 cases
                                      of the iteration with the lowest CAD/WCAD
                                      score.
      'batch_size':[x] - no. of iterations in a batch
      'plncases':[x,x,...] - the daily number of local COVID-19 cases
                             from empical data    lowest CAD/WCAD score.
      'pdates':[ddd,ddd,...] - the date & time corresponding to the days in
                              'inf_days_bins'.
      'seed':[xxxx] - The seed used to prime NumPy's SeedSequence.
    '''
    data = {}
    with open(source, mode='r', encoding='utf-8') as csvfile:
        csv_reader = reader(csvfile, delimiter=',')
        header = next(csv_reader)
        for n, (k, v) in enumerate(zip(header, csv_reader)):
            if n != 8:
                if n != 1:
                    v = [int(i) for i in v]
                else:
                    if v[1].isnumeric() is False:
                        v = [int(v[0]), float(v[1])]
                    else:
                        v = [int(i) for i in v]
            else:
                v = [convert_str_to_datetime(
                    i, dateformat='%Y-%m-%d %H:%M:%S%z') for i in v]
            data[k] = v
    return data


def plot_CAD_WCAD(sub_dirs, mus_range, cad_scores, wcad_scores):
    fontsize_axis = 16
    fontsize_labels = 20
    markerstyles = ['d', 'x', 'o']
    labels = ['Seed1', 'Seed2', 'Seed3']
    colors = ['r', 'g', 'b']

    fig, ax = plt.subplots(1, 2, sharex=True, sharey=False)
    for n, dir in enumerate(sub_dirs):
        ax[0].plot(
            mus_range, cad_scores[dir.name], label=labels[n], alpha=0.8,
            lw=4, marker=markerstyles[n], ms=10, color=colors[n])
        ax[1].plot(
            mus_range, wcad_scores[dir.name], label=labels[n], alpha=0.8,
            lw=4, marker=markerstyles[n], ms=10, color=colors[n])
    # Configure Y-axis
    ax[0].set_ylabel('CAD', fontweight='bold', fontsize=fontsize_labels,
                     fontstyle='italic')
    ax[1].set_ylabel('WCAD', fontweight='bold', fontsize=fontsize_labels,
                     fontstyle='italic')
    ax[0].tick_params(axis='both', which='both', labelsize=fontsize_axis)
    ax[1].tick_params(axis='both', which='both', labelsize=fontsize_axis)
    # Configure X-axis
    xlabel = r'Step0 Unselected Missing $\mu_{mean}$ elements values'
    ax[0].set_xlabel(xlabel, fontweight='bold', fontsize=fontsize_labels)
    ax[1].set_xlabel(xlabel, fontweight='bold', fontsize=fontsize_labels)
    ax[0].xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax[0].set_xlim([min(mus_range)-1, max(mus_range)+1])
    xticklabels = [str(i) for i in range(len(mus_range)+2)]
    xticklabels[0] = xticklabels[-1] = ''
    ax[0].set_xticks(np.arange(len(mus_range)+2), labels=xticklabels)
    # Draw Grid
    ax[0].grid(linestyle=':', color='grey')
    ax[1].grid(linestyle=':', color='grey')
    # Draw Legend
    ax[0].legend(loc='upper right', fontsize=fontsize_labels)
    ax[1].legend(loc='upper right', fontsize=fontsize_labels)
    # Configure Subplot layput
    top = 0.950  # the top of the subplots of the figure
    bottom = 0.090  # the bottom of the subplots of the figure
    left = 0.070  # the left side of the subplots of the figure
    right = 0.970  # the right side of the subplots of the figure
    hspace = 0.163  # the amount of height reserved for space between subplots,
    # expressed as a fraction of the average axis height
    wspace = 0.175  # the amount of width reserved for space between subplots,
    # expressed as a fraction of the average axis width
    fig.subplots_adjust(left=left, bottom=bottom, right=right, top=top,
                        wspace=wspace, hspace=hspace)
    plt.show()


def get_min_cad_wcad_from_1step(sub_dirs, mus_range):
    # Define dictionary to store the CAD and WCAD results
    # CADS = {'s1':[r1,r2,....], 's2':[r1,r2,....], 's3':[r1,r2,....]}
    # WCADS = {'s1':[r1,r2,....], 's2':[r1,r2,....], 's3':[r1,r2,....]}
    # r1, r2, denotes the results
    CADS = {}
    WCADS = {}
    cad_scores = {}
    wcad_scores = {}
    min_cad_scores = {}
    min_wcad_scores = {}
    for nR, dir in enumerate(sub_dirs):
        # Store CAD and WCAD data
        CADS[dir.name] = {mu: read_CAD_WCAD(
            dir / f'{dir.name}_stp0_o{mu}_cad.csv')
            for mu in mus_range
        }
        WCADS[dir.name] = {mu: read_CAD_WCAD(
            dir / f'{dir.name}_stp0_o{mu}_wcad.csv')
            for mu in mus_range
        }
        # Store CAD and WCAD scores
        cad_scores[dir.name] = np.array(
            [CADS[dir.name][i]['score'][1] for i in mus_range])
        wcad_scores[dir.name] = np.array(
            [WCADS[dir.name][i]['score'][1] for i in mus_range])
        # Store Index of lowest CAD and WCAD scores
        scores = cad_scores[dir.name]
        min_cad_scores[dir.name] = (
            np.amin(scores), np.where(scores == np.amin(scores))[0].item())
        scores = wcad_scores[dir.name]
        min_wcad_scores[dir.name] = (
            np.amin(scores), np.where(scores == np.amin(scores))[0].item())
    for dir in sub_dirs:
        for i in mus_range:
            print(dir.name, i, CADS[dir.name][i])

    for dir in sub_dirs:
        for i in mus_range:
            print(dir.name, i, WCADS[dir.name][i])
    print()
    for dir in sub_dirs:
        print(dir.name, cad_scores[dir.name])
    print()
    for dir in sub_dirs:
        print(dir, wcad_scores[dir.name])
    print()
    print(min_cad_scores)
    print(min_wcad_scores)
    return CADS, WCADS, cad_scores, wcad_scores, min_cad_scores,\
        min_wcad_scores


if __name__ == "__main__":
    # Set up results directories path
    repo_dir = Path(__file__).resolve().parents[2]
    results_data = repo_dir / '2_Results' / '2_μmean_COVID19_SARSCOV2_trends'
    s1 = results_data / 's1'
    s2 = results_data / 's2'
    s3 = results_data / 's3'
    sub_dirs = [s1, s2, s3]
    # Plot results
    mus_range = [i for i in range(1, 19, 1)]
    cads, wcads, cad_scores, wcad_scores, min_cad_scores, min_wcad_scores = \
        get_min_cad_wcad_from_1step(sub_dirs, mus_range)
    # Plot CAD & WCAD results
    plot_CAD_WCAD(sub_dirs, mus_range, cad_scores, wcad_scores)
