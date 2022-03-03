#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import MO, TU, WE, TH, FR, SA, SU
from matplotlib.patches import Rectangle
from pathlib import Path
from datetime import datetime, timedelta, timezone
import csv
from debug_functions import prv, pri
import numpy as np
import statistics as stats


def read_CAD_WCAD(source='source.csv'):
    '''Method to extract data from either the CAD.csv and WCAD.csv result file
    of the Predict_SARSCOV2_COVID19_VIA_CPTPP class and return these data as a
    dictionary.

    Dictionary keys and values (note: each value  in values is a list obj):
      'batch':[xxx] - give the batch number (integer) that gave the lowest CAD
                      or WCAD score.
      'least_score':[iii, iii] - gives the iteration number of that batch and
                                 the CAD/WCAD score from that iteration.
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
      'plncases':[x,x,...] - the daily number of local COVID-19 cases from
                             empirical data of the lowest CAD/WCAD score.
      'pdates':[ddd,ddd,...] - the date & time corresponding to the days in
                               'inf_days_bins'.
      'seed':[xxxx] - The seed used to prime NumPy's SeedSequence.
    '''
    data = {}
    with open(source, mode='r', encoding='utf-8') as csvfile:
        csv_reader = csv.reader(csvfile, delimiter=',')
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


def convert_str_to_datetime(date, dateformat="%d/%m/%Y", tzinfo=timezone.utc):
    refdate = datetime.strptime(date, dateformat)
    refdate = refdate.replace(tzinfo=tzinfo)
    return refdate


def get_data(dictionary):
    cba_pdates = []
    cbs_pdates = []
    cbp1_pdates = []
    cbp2_pdates = []
    cba_cy = []
    cbs_cy = []
    cbp1_cy = []
    cbp2_cy = []
    emp_cba_cy = []
    emp_cbs_cy = []
    emp_cbp1_cy = []
    emp_cbp2_cy = []
    for n, (seed, v) in enumerate(dictionary.items()):
        cba = convert_str_to_datetime('03/04/2020')  # cb announced
        cbs = convert_str_to_datetime('07/04/2020')  # cb starts
        cbp1 = convert_str_to_datetime('20/04/2020')  # cb peak1
        cbp2 = convert_str_to_datetime('01/05/2020')  # cb peak2
        pdates = v['pdates']
        emp = np.array(v['plncases'], dtype=np.int16)
        sarscov2 = np.array(v['inf_days_ncases'], dtype=np.int16)
        cba_index = pdates.index(cba)
        cbs_index = pdates.index(cbs)
        cbp1_index = pdates.index(cbp1)
        cbp2_index = pdates.index(cbp2)
        cba_pdates.append(pdates[:cba_index+1])
        cbs_pdates.append(pdates[:cbs_index+1])
        cbp1_pdates.append(pdates[:cbp1_index+1])
        cbp2_pdates.append(pdates[:cbp2_index+1])
        cba_cy.append(sarscov2[:cba_index+1].cumsum())
        cbs_cy.append(sarscov2[:cbs_index+1].cumsum())
        cbp1_cy.append(sarscov2[:cbp1_index+1].cumsum())
        cbp2_cy.append(sarscov2[:cbp2_index+1].cumsum())
        emp_cba_cy.append(emp[:cba_index+1].cumsum())
        emp_cbs_cy.append(emp[:cbs_index+1].cumsum())
        emp_cbp1_cy.append(emp[:cbp1_index+1].cumsum())
        emp_cbp2_cy.append(emp[:cbp2_index+1].cumsum())
    return (cba_pdates, cbs_pdates, cbp1_pdates, cbp2_pdates, cba_cy, cbs_cy,
            cbp1_cy, cbp2_cy, emp_cba_cy, emp_cbs_cy, emp_cbp1_cy, emp_cbp2_cy)


def plot_cumulativeSARSCoV2(seeds, cads_data, wcads_data):
    # Plot Figure
    fontsize_axis = 16
    fontsize_axislabel = 20
    fig, ax = plt.subplots(2, 2, sharex=True, sharey=False)
    axes = (ax[0, 0], ax[1, 0], ax[0, 1], ax[1, 1])
    colors = ('r', 'g', 'b', 'C4', 'C5', 'C6')
    marks = ('d', 'x', 'o')

    def add_figure_labels():
        msgs = ['(a) At CB Announcement', '(b) At Start Of CB',
                r'(c) At $1^{st}$ SARS-CoV-2 Peak', r'(d) At $2^{nd}$ SARS-CoV-2 Peak']
        msgx = mdates.date2num(convert_str_to_datetime('26/01/2020'))
        msgys = [250, 500, 1000, 2000]
        for n, msg in enumerate(msgs):
            axes[n].text(msgx, msgys[n], msg, ha="left", va="bottom",
                         color='black', fontsize=fontsize_axislabel, alpha=0.8)

    def add_circuit_breaker(ax, subplot, ymin, ymax):
        cb_ann_x = mdates.date2num(convert_str_to_datetime('03/04/2020'))
        ax.vlines(cb_ann_x, 0, 1.0, transform=ax.get_xaxis_transform(),
                  colors='orange', linestyles='solid', alpha=0.8)
        cb_start = mdates.date2num(convert_str_to_datetime('07/04/2020'))
        cb_end = mdates.date2num(convert_str_to_datetime('01/06/2020'))
        width = cb_end - cb_start
        # Boundary
        height = (ymax + abs(ymin)) * 1.2
        rect = Rectangle((cb_start, ymin), width, height, linewidth=0,
                         edgecolor=None, facecolor='grey', alpha=0.4)
        ax.add_patch(rect)
        if subplot == 1:
            msg = 'Circuit Breaker (CB)'
            msgx = mdates.date2num(convert_str_to_datetime('26/04/2020'))
            msgy = 0
            ax.text(msgx, msgy, msg, ha="left", va="bottom", color='black',
                    fontsize=fontsize_axislabel, alpha=0.8, rotation=90)
        if subplot == 3:
            msg = r'Circuit Breaker Announcement.'
            msgx = mdates.date2num(convert_str_to_datetime('01/04/2020'))
            msgy = 5000
            ax.text(msgx, msgy, msg, ha="left", va="bottom", color="blue",
                    fontsize=12, rotation=90, alpha=0.8)

    cads_sars_xs = cads_data[0:4]
    cads_sars_ys = cads_data[4:8]
    cads_emp_ys = cads_data[8:12]
    wcads_sars_xs = wcads_data[0:4]
    wcads_sars_ys = wcads_data[4:8]
    wcads_emp_ys = wcads_data[8:12]

    # for i in cads_sars_xs:
    #     print(type(i), len(i), i)
    #     for j in i:
    #         print(type(j), len(j), j)

    # cba_pdates, cbs_pdates, cbp1_pdates, cbp2_pdates, cba_cy, cbs_cy,
    # cbp1_cy, cbp2_cy, emp_cba_cy, emp_cbs_cy, emp_cbp1_cy, emp_cbp2_cy
    # CAD
    parliment_by_lw = convert_str_to_datetime('25/03/2020')
    parliment_SARSCOV2_population = []
    parliment_COVID19_population = []
    label = 'SARS-CoV-2:CAD-'
    for n, xdates in enumerate(cads_sars_xs[0]):
        axes[0].plot(xdates, cads_sars_ys[0][n], label=label+f'{seeds[n]}',
                     marker=marks[n], color=colors[n], alpha=0.8, lw=0)
        parliment_index = xdates.index(parliment_by_lw)
        spop = cads_sars_ys[0][n][parliment_index]
        parliment_SARSCOV2_population.append(spop)
        print('CAD', parliment_by_lw, spop)
    for n, xdates in enumerate(cads_sars_xs[1]):
        axes[1].plot(xdates, cads_sars_ys[1][n], label=label+f'{seeds[n]}',
                     marker=marks[n], color=colors[n], alpha=0.8, lw=0)
    for n, xdates in enumerate(cads_sars_xs[2]):
        axes[2].plot(xdates, cads_sars_ys[2][n], label=label+f'{seeds[n]}',
                     marker=marks[n], color=colors[n], alpha=0.8, lw=0)
    for n, xdates in enumerate(cads_sars_xs[3]):
        axes[3].plot(xdates, cads_sars_ys[3][n], label=label+f'{seeds[n]}',
                     marker=marks[n], color=colors[n], alpha=0.8, lw=0)
    # WCAD
    label = 'SARS-CoV-2:WCAD-'
    for n, xdates in enumerate(wcads_sars_xs[0]):
        axes[0].plot(xdates, wcads_sars_ys[0][n], label=label+f'{seeds[n]}',
                     marker=marks[n], color=colors[n+3], alpha=0.8, lw=0)
        parliment_index = xdates.index(parliment_by_lw)
        spop = cads_sars_ys[0][n][parliment_index]
        parliment_SARSCOV2_population.append(spop)
        print('WCAD', parliment_by_lw, spop)
    for n, xdates in enumerate(wcads_sars_xs[1]):
        axes[1].plot(xdates, wcads_sars_ys[1][n], label=label+f'{seeds[n]}',
                     marker=marks[n], color=colors[n+3], alpha=0.8, lw=0)
    for n, xdates in enumerate(wcads_sars_xs[2]):
        axes[2].plot(xdates, wcads_sars_ys[2][n], label=label+f'{seeds[n]}',
                     marker=marks[n], color=colors[n+3], alpha=0.8, lw=0)
    for n, xdates in enumerate(wcads_sars_xs[3]):
        axes[3].plot(xdates, wcads_sars_ys[3][n], label=label+f'{seeds[n]}',
                     marker=marks[n], color=colors[n+3], alpha=0.8, lw=0)
    # Empirical Covid-19
    label = 'COVID-19:Empirical'
    axes[0].plot(cads_sars_xs[0][0], cads_emp_ys[0][0], label=label,
                 marker='s', color='black', alpha=0.8, lw=2)
    parliment_index = cads_sars_xs[0][0].index(parliment_by_lw)
    cpop = cads_emp_ys[0][0][parliment_index]
    parliment_COVID19_population.append(cpop)
    print('COVID19', parliment_by_lw, cpop)
    axes[1].plot(cads_sars_xs[1][0], cads_emp_ys[1][0], label=label,
                 marker='s', color='black', alpha=0.8, lw=2)
    axes[2].plot(cads_sars_xs[2][0], cads_emp_ys[2][0], label=label,
                 marker='s', color='black', alpha=0.8, lw=2)
    axes[3].plot(cads_sars_xs[3][0], cads_emp_ys[3][0], label=label,
                 marker='s', color='black', alpha=0.8, lw=2)
    # Analyse Parliment day Population
    pri('\nparliment_SARSCOV2_population', parliment_SARSCOV2_population)
    print(stats.mean(parliment_SARSCOV2_population), stats.stdev(parliment_SARSCOV2_population))
    pri('parliment_COVID19_population', parliment_COVID19_population)
    # Configure Subplots 1 & 2
    ymax = [max(cads_sars_ys[0][0].max(), cads_sars_ys[0][1].max(), cads_sars_ys[0][2].max()),
            max(cads_sars_ys[1][0].max(), cads_sars_ys[1][1].max(), cads_sars_ys[1][2].max()),
            max(cads_sars_ys[2][0].max(), cads_sars_ys[2][1].max(), cads_sars_ys[2][2].max()),
            max(cads_sars_ys[3][0].max(), cads_sars_ys[3][1].max(), cads_sars_ys[3][2].max())]
    ymin = [-50, -100, -750, -1000]
    for n, ax in enumerate(axes):
        # Add circuit-breaker zone
        add_circuit_breaker(ax, n, ymin[n], ymax[n])
        # Configure Y-axis
        ylabel = r'Population'
        ax.set_ylabel(ylabel, fontweight='bold', fontsize=fontsize_axislabel)
        ax.tick_params(axis='both', which='both', labelsize=fontsize_axis)
        # ax.yaxis.set_major_locator(ticker.MultipleLocator(2))
        # ax.yaxis.set_minor_locator(ticker.MultipleLocator(1))
        ax.set_ylim([ymin[n], (abs(ymin[n])+ymax[n])*1.1])
        # Draw Grid
        ax.grid(linestyle=':', color='grey')
        # Draw Legend
        ax.legend(loc='upper left', fontsize=fontsize_axislabel)
        # Configure X-axis
        ax.get_xaxis().set_label_text(label='Press Release Dates',
                                      fontsize=fontsize_axislabel,
                                      fontweight='bold')
        fig.autofmt_xdate()  # rotate and align the tick labels so they look better. From https://matplotlib.org/3.1.3/gallery/recipes/common_date_problems.html
        ax.xaxis.set_major_locator(
            mdates.WeekdayLocator(byweekday=SU, interval=1, tz=None))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d %a", tz=None))
        ax.xaxis.set_minor_locator(
            mdates.WeekdayLocator(byweekday=(MO, TU, WE, TH, FR, SA),
                                  interval=1, tz=None))
        # ax.xaxis.set_minor_formatter(mdates.DateFormatter("%a", tz=None))
        # ax.set_xlim(auto=True)
        start = convert_str_to_datetime('23/01/2020')
        start = start - timedelta(days=2.)
        end = convert_str_to_datetime('01/05/2020')
        end = end + timedelta(days=2.)
        ax.set_xlim([start, end])
    # Add Subplot labels
    add_figure_labels()
    # Rotate xaxis minorticklabels
    # manual control
    angle = 90
    plt.setp(axes[1].get_xmajorticklabels(), rotation=angle, ha='center')
    plt.setp(axes[3].get_xmajorticklabels(), rotation=angle, ha='center')
    # Configure Subplot layput
    top = 0.980  # the top of the subplots of the figure
    bottom = 0.130  # the bottom of the subplots of the figure
    left = 0.060  # the left side of the subplots of the figure
    right = 0.985  # the right side of the subplots of the figure
    hspace = 0.050  # the amount of height reserved for space between subplots,
    # expressed as a fraction of the average axis height
    wspace = 0.220  # the amount of width reserved for space between subplots,
    # expressed as a fraction of the average axis width
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top,
                        wspace=wspace, hspace=hspace)
    cymax = [cads_sars_ys[0][0].max(), cads_sars_ys[0][1].max(), cads_sars_ys[0][2].max(),
             cads_sars_ys[1][0].max(), cads_sars_ys[1][1].max(), cads_sars_ys[1][2].max(),
             cads_sars_ys[2][0].max(), cads_sars_ys[2][1].max(), cads_sars_ys[2][2].max(),
             cads_sars_ys[3][0].max(), cads_sars_ys[3][1].max(), cads_sars_ys[3][2].max()]
    wymax = [wcads_sars_ys[0][0].max(), wcads_sars_ys[0][1].max(), wcads_sars_ys[0][2].max(),
             wcads_sars_ys[1][0].max(), wcads_sars_ys[1][1].max(), wcads_sars_ys[1][2].max(),
             wcads_sars_ys[2][0].max(), wcads_sars_ys[2][1].max(), wcads_sars_ys[2][2].max(),
             wcads_sars_ys[3][0].max(), wcads_sars_ys[3][1].max(), wcads_sars_ys[3][2].max()]
    covidmax = [cads_emp_ys[0][0].max(), cads_emp_ys[1][0].max(),
                cads_emp_ys[2][0].max(), cads_emp_ys[3][0].max()]
    pri('cymax', cymax)
    pri('wymax', wymax)
    pri('covidmax', covidmax)
    plt.show()


def get_cumlulativedifference(dictionary):
    cme = []
    sme = []
    smc = []
    for n, (seed, v) in enumerate(dictionary.items()):
        pdates = v['pdates']
        emp = np.array(v['plncases'], dtype=np.int16).cumsum()
        covid19 = np.array(v['covid_days_ncases'], dtype=np.int16).cumsum()
        sarscov2 = np.array(v['inf_days_ncases'], dtype=np.int16).cumsum()
        cme.append(covid19 - emp)
        sme.append(sarscov2 - emp)
        smc.append(sarscov2 - covid19)
    for i in cme:
        pri('cme', i)
        print('total = ', i.sum())
    for i in sme:
        pri('sme', i)
        print('total = ', i.sum())
    for i in smc:
        pri('smc', i)
        print('total = ', i.sum())
    return pdates, cme, sme, smc


def plot_cumlulativedifference(seeds, cads_data, wcads_data):
    # Plot Figure
    fontsize_axis = 16
    fontsize_axislabel = 20
    fig, ax = plt.subplots(3, 1, sharex=True, sharey=False)
    axes = (ax[0], ax[1], ax[2])
    colors = ('r', 'g', 'b', 'C4', 'C5', 'C6')
    marks = ('d', 'x', 'o', '^', '+', 's')

    def add_circuit_breaker(ax, subplot, ymin, ymax):
        cb_ann_x = mdates.date2num(convert_str_to_datetime('03/04/2020'))
        ax.vlines(cb_ann_x, 0, 1.0, transform=ax.get_xaxis_transform(),
                  colors='orange', linestyles='solid', alpha=0.5)
        cb_start = mdates.date2num(convert_str_to_datetime('07/04/2020'))
        cb_end = mdates.date2num(convert_str_to_datetime('01/06/2020'))
        width = cb_end - cb_start
        # Boundary
        height = (ymax + abs(ymin)) * 2.
        rect = Rectangle((cb_start, ymin), width, height, linewidth=0,
                         edgecolor=None, facecolor='grey', alpha=0.2)
        ax.add_patch(rect)
        if subplot == 1:
            msg = 'Circuit Breaker (CB)'
            msgx = mdates.date2num(convert_str_to_datetime('26/04/2020'))
            msgy = 2000
            ax.text(msgx, msgy, msg, ha="left", va="bottom", color='black',
                    fontsize=fontsize_axislabel, alpha=0.5, rotation=0)
            msg = r'Circuit Breaker Announcement.'
            msgx = mdates.date2num(convert_str_to_datetime('01/04/2020'))
            msgy = 0
            ax.text(msgx, msgy, msg, ha="left", va="bottom", color="blue",
                    fontsize=12, rotation=90, alpha=0.6)
    # CAD
    for n, data in enumerate(cads_data[1]):
        axes[0].plot(cads_data[0], data, label=f'CAD-{seeds[n]}',
                     marker=marks[n], color=colors[n], alpha=0.3, lw=0)
    for n, data in enumerate(cads_data[2]):
        axes[1].plot(cads_data[0], data, label=f'CAD-{seeds[n]}',
                     marker=marks[n], color=colors[n], alpha=0.3, lw=0)
    for n, data in enumerate(cads_data[3]):
        axes[2].plot(cads_data[0], data, label=f'CAD-{seeds[n]}',
                     marker=marks[n], color=colors[n], alpha=0.3, lw=0)
    # WCADS
    for n, data in enumerate(wcads_data[1]):
        axes[0].plot(wcads_data[0], data, label=f'WCAD-{seeds[n]}',
                     marker=marks[n+3], color=colors[n+3], alpha=0.3, lw=0)
    for n, data in enumerate(wcads_data[2]):
        axes[1].plot(wcads_data[0], data, label=f'WCAD-{seeds[n]}',
                     marker=marks[n+3], color=colors[n+3], alpha=0.3, lw=0)
    for n, data in enumerate(wcads_data[3]):
        axes[2].plot(wcads_data[0], data, label=f'WCAD-{seeds[n]}',
                     marker=marks[n+3], color=colors[n+3], alpha=0.3, lw=0)
    # Configure Subplots
    ymax = [500, 5000, 5000]
    ymin = [-750, -100, -100]
    for n, ax in enumerate(axes):
        # Add circuit-breaker zone
        add_circuit_breaker(ax, n, ymin[n], ymax[n])
        # Configure Y-axis
        ylabel = r'Number of Cases'
        ax.set_ylabel(ylabel, fontweight='bold', fontsize=fontsize_axislabel)
        ax.tick_params(axis='both', which='both', labelsize=fontsize_axis)
        # ax.yaxis.set_major_locator(ticker.MultipleLocator(2))
        # ax.yaxis.set_minor_locator(ticker.MultippopulationleLocator(1))
        ax.set_ylim([ymin[n], (abs(ymin[n])+ymax[n])*1.1])
        # Draw Grid
        ax.grid(linestyle=':', color='grey')
        # Draw Legend
        ax.legend(loc='upper left', fontsize=fontsize_axislabel)
        # Configure X-axis
        ax.get_xaxis().set_label_text(label='Press Release Dates',
                                      fontsize=fontsize_axislabel,
                                      fontweight='bold')
        fig.autofmt_xdate()  # rotate and align the tick labels so they look better. From https://matplotlib.org/3.1.3/gallery/recipes/common_date_problems.html
        ax.xaxis.set_major_locator(
            mdates.WeekdayLocator(byweekday=SU, interval=1, tz=None))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d %a", tz=None))
        ax.xaxis.set_minor_locator(
            mdates.WeekdayLocator(byweekday=(MO, TU, WE, TH, FR, SA),
                                  interval=1, tz=None))
        # ax.xaxis.set_minor_formatter(mdates.DateFormatter("%a", tz=None))
        # ax.set_xlim(auto=True)
        start = cads_data[0][0]
        start = start - timedelta(days=2.)
        end = cads_data[0][-1]
        end = end + timedelta(days=2.)
        ax.set_xlim([start, end])
    axes[0].set_ylabel(r'$P_{COVID-19}-E_{COVID-19}$', fontweight='bold',
                       fontsize=fontsize_axislabel)
    axes[1].set_ylabel(r'$P_{SARS-CoV-2}-E_{COVID-19}$', fontweight='bold',
                       fontsize=fontsize_axislabel)
    axes[2].set_ylabel(r'$P_{SARS-CoV-2}-P_{COVID-19}$', fontweight='bold',
                       fontsize=fontsize_axislabel)
    # Rotate xaxis minorticklabels
    # manual control
    angle = 90
    # plt.setp(axes[0].get_xmajorticklabels(), rotation=angle, ha='center')
    plt.setp(axes[2].get_xmajorticklabels(), rotation=angle, ha='center')
    # Configure Subplot layput
    top = 0.980  # the top of the subplots of the figure
    bottom = 0.130  # the bottom of the subplots of the figure
    left = 0.060  # the left side of the subplots of the figure
    right = 0.985  # the right side of the subplots of the figure
    hspace = 0.050  # the amount of height reserved for space between subplots,
    # expressed as a fraction of the average axis height
    wspace = 0.220  # the amount of width reserved for space between subplots,
    # expressed as a fraction of the average axis width
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top,
                        wspace=wspace, hspace=hspace)


if __name__ == "__main__":
    # Results from Modified Resemblance Alogorithm Analyses
    repo_dir = Path(__file__).resolve().parents[2]
    results_data = repo_dir / '2_Results' / '2_μmean_COVID19_SARSCOV2_trends'
    source1 = results_data / 's1' / 's1_stp2_cad_cad.csv'
    source2 = results_data / 's2' / 's2_stp2_cad_cad.csv'
    source3 = results_data / 's3' / 's3_stp2_cad_cad.csv'
    cadfiles = (source1, source2, source3)
    source4 = results_data / 's1' / 's1_stp2_wcad_wcad.csv'
    source5 = results_data / 's2' / 's2_stp2_wcad_wcad.csv'
    source6 = results_data / 's3' / 's3_stp2_wcad_wcad.csv'
    wcadfiles = (source4, source5, source6)
    seeds = ('Seed1', 'Seed2', 'Seed3')
    CADS = {seeds[n]: read_CAD_WCAD(source=i) for n, i in enumerate(cadfiles)}
    WCADS = {seeds[n]: read_CAD_WCAD(source=i) for n, i in enumerate(wcadfiles)}
    cads_data = get_data(CADS)
    wcads_data = get_data(WCADS)
    plot_cumulativeSARSCoV2(seeds, cads_data, wcads_data)
    # cads_data = get_cumlulativedifference(CADS)
    # wcads_data = get_cumlulativedifference(WCADS)
    # plot_cumlulativedifference(seeds, cads_data, wcads_data)
    plt.show()
