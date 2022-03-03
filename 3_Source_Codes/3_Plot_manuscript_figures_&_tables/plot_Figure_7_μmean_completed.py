import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
from matplotlib.dates import MO, TU, WE, TH, FR, SA, SU
from matplotlib.patches import Rectangle
from pathlib import Path
from datetime import datetime, timedelta, timezone
import csv
from debug_functions import prv, pri
import numpy as np


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


def get_sma(xdate, y, window):
    # Get the Moving Average of predicted Mus ( 14 days window used )
    # xsma = np.arange( max(xmean) )
    def get_moving_average(values):
        weights = np.repeat(1.0, window) / window
        smas = np.convolve(values, weights, 'valid')
        return smas.tolist()
    M = window
    N = len(y)
    size = max(M, N) - min(M, N) + 1
    prv('size', size)
    edge = window // 2
    prv('edge', edge)
    if window % 2 == 0:
        # even
        print(' even')
        xdate_sma = xdate[edge:-(edge-1)]
    else:
        # odd
        print(' odd')
        xdate_sma = xdate[edge:-(edge)]
    y_sma = get_moving_average(y)
    pri('xdate_sma', xdate_sma)
    pri('y_sma', y_sma)
    return xdate_sma, y_sma


def plot_pmus_curves(CADS, WCADS):
    fontsize_axis = 16
    fontsize_axislabel = 20

    def add_circuit_breaker(ax, subplot):
        lcovid0_x = mdates.date2num(convert_str_to_datetime('02/02/2020'))
        ax.vlines(lcovid0_x, 0, 1.0, transform=ax.get_xaxis_transform(),
                  colors='orange', linestyles='solid', alpha=0.8)
        cb_ann_x = mdates.date2num(convert_str_to_datetime('03/04/2020'))
        ax.vlines(cb_ann_x, 0, 1.0, transform=ax.get_xaxis_transform(),
                  colors='orange', linestyles='solid', alpha=0.8)
        cb_start = mdates.date2num(convert_str_to_datetime('07/04/2020'))
        cb_end = mdates.date2num(convert_str_to_datetime('01/06/2020'))
        width = cb_end - cb_start
        # Boundary
        rect = Rectangle((cb_start, 0), width, 19, linewidth=0,
                         edgecolor=None, facecolor='grey', alpha=0.4)
        ax.add_patch(rect)
        if subplot == 1:
            # Text
            msg = 'Circuit Breaker'
            fontcolor = 'black'
            msgx = mdates.date2num(convert_str_to_datetime('03/05/2020'))
            msgy = 16
            ax.text(msgx, msgy, msg, ha="left", va="bottom", color=fontcolor,
                    fontsize=fontsize_axislabel, alpha=0.8)
            msgx = mdates.date2num(convert_str_to_datetime('01/04/2020'))
            msg = r'$3^{rd}$ Apr: Circuit Breaker Announced.'
            ax.text(msgx, 10, msg, ha="left", va="bottom", color="blue",
                    fontsize=12, rotation=90, alpha=0.8)
        if subplot == 0:
            # Text
            msgx = mdates.date2num(convert_str_to_datetime('31/01/2020'))
            msg = r'$2^{nd}$ Feb: $1^{st}$ Local COVID-19.'
            ax.text(msgx, 0.2, msg, ha="left", va="bottom", color="blue",
                    fontsize=12, rotation=90, alpha=0.8)

    # Plot Figure
    fig, ax = plt.subplots(2, 1, sharex=True, sharey=False)
    axes = (ax[0], ax[1])
    colors = ('r', 'g', 'b')
    marks = ('d', 'x', 'o')
    # Plot subplots 1
    for n, (seed, v) in enumerate(CADS.items()):
        x_dates = [mdates.date2num(date) for date in v['pdates']]
        y_pmus = v['final_mus']
        score = v['score'][1]
        # Plot data
        axes[0].plot(x_dates, y_pmus, label=f'{seed}, CAD={score}',
                     marker=marks[n], color=colors[n], alpha=0.8, lw=2)
        window = 14
        fxsma, fysma = get_sma(x_dates, y_pmus, window)
        ax[0].plot(fxsma, fysma, label=f'{seed}-Moving Average(14-days)',
                   color='C7', alpha=0.8, lw=4)
    # Plot subplots 2
    for n, (seed, v) in enumerate(WCADS.items()):
        x_dates = [mdates.date2num(date) for date in v['pdates']]
        y_pmus = v['final_mus']
        score = round(float(v['score'][1]), 2)
        # Plot data
        axes[1].plot(x_dates, y_pmus, label=f'{seed}, WCAD={score}',
                     marker=marks[n], color=colors[n], alpha=0.8, lw=2)
        window = 14
        fxsma, fysma = get_sma(x_dates, y_pmus, window)
        ax[1].plot(fxsma, fysma, label=f'{seed}-Moving Average(14-days)',
                   color='C7', alpha=0.8, lw=4)
    # Configure Subplots 1 & 2
    for n in range(2):
        # Configure Y-axis
        ylabel = r'Mean of the $\mu$ Estimates, $\mu_{mean}$'
        axes[n].set_ylabel(ylabel, fontweight='bold',
                           fontsize=fontsize_axislabel)
        axes[n].tick_params(axis='both', which='both', labelsize=fontsize_axis)
        axes[n].yaxis.set_major_locator(ticker.MultipleLocator(2))
        axes[n].yaxis.set_minor_locator(ticker.MultipleLocator(1))
        axes[n].set_ylim([0, 19])
        # Draw Grid
        ax[n].grid(linestyle=':', color='grey')
        # Draw Legend
        ax[n].legend(loc='upper left', fontsize=fontsize_axislabel)
        # Draw Circuit Breaker
        add_circuit_breaker(ax[n], n)
    # Configure X-axis
    ax[1].get_xaxis().set_label_text(label='Dates',
                                     fontsize=fontsize_axislabel,
                                     fontweight='bold')
    fig.autofmt_xdate()  # rotate and align the tick labels so they look better. From https://matplotlib.org/3.1.3/gallery/recipes/common_date_problems.html
    ax[1].xaxis.set_major_locator(
        mdates.WeekdayLocator(byweekday=SU, interval=1, tz=None))
    ax[1].xaxis.set_major_formatter(mdates.DateFormatter("%b %d %a", tz=None))
    ax[1].xaxis.set_minor_locator(
        mdates.WeekdayLocator(byweekday=(MO, TU, WE, TH, FR, SA),
                              interval=1, tz=None))
    # ax.xaxis.set_minor_formatter(mdates.DateFormatter("%a", tz=None))
    # ax.set_xlim(auto=True)
    start = convert_str_to_datetime('23/01/2020')
    start = start - timedelta(days=2.)
    end = convert_str_to_datetime('21/08/2020')
    end = end + timedelta(days=2.)
    ax[1].set_xlim([start, end])
    # Rotate xaxis minorticklabels
    # manual control
    angle = 90
    plt.setp(ax[1].get_xminorticklabels(), rotation=angle, ha='center')
    plt.setp(ax[1].get_xmajorticklabels(), rotation=angle, ha='center')
    handles0, labels0 = ax[0].get_legend_handles_labels()
    handles1, labels1 = ax[1].get_legend_handles_labels()
    print(handles0)
    print(labels0)
    handles0 = [handles0[0], handles0[2], handles0[4], handles0[1]]
    handles1 = [handles1[0], handles1[2], handles1[4], handles1[1]]
    labels0 = [labels0[0], labels0[2], labels0[4], 'Moving Average(14-days)']
    labels1 = [labels1[0], labels1[2], labels1[4], 'Moving Average(14-days)']
    ax[0].legend(handles0, labels0, loc='upper left', fontsize=fontsize_axislabel)
    ax[1].legend(handles1, labels1, loc='upper left', fontsize=fontsize_axislabel)
    # Configure Subplot layput
    top = 0.980  # the top of the subplots of the figure
    bottom = 0.130  # the bottom of the subplots of the figure
    left = 0.050  # the left side of the subplots of the figure
    right = 0.985  # the right side of the subplots of the figure
    hspace = 0.050  # the amount of height reserved for space between subplots,
    # expressed as a fraction of the average axis height
    wspace = 0.220  # the amount of width reserved for space between subplots,
    # expressed as a fraction of the average axis width
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top,
                        wspace=wspace, hspace=hspace)
    plt.show()


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
    plot_pmus_curves(CADS, WCADS)
