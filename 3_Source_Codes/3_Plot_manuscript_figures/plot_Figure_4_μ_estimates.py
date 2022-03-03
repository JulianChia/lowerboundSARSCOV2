import numpy as np

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
from matplotlib.dates import MO, TU, WE, TH, FR, SA, SU
from matplotlib.patches import Rectangle

from pathlib import Path
from datetime import datetime, timedelta, timezone
import statistics as stats
from math import ceil
import csv

fontsize_axis = 16
fontsize_axislabel = 20


def read_predicted_mus(source='out.csv'):
    '''Method to extract predicted mus from a .csv file.'''
    predictions = {}
    dates = []
    with open(source, mode='r', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        next(reader)  # Skip header row
        for row in reader:
            if row:
                day = int(row[0])
                date = convert_str_to_datetime(row[1], dateformat='%Y-%m-%d')
                predictions[day] = [int(i) for i in row[2:]]
                dates.append(date)
    return predictions, dates


def convert_str_to_datetime(date, dateformat="%d/%m/%Y", tzinfo=timezone.utc):
    refdate = datetime.strptime(date, dateformat)
    refdate = refdate.replace(tzinfo=tzinfo)
    return refdate


def get_mean(y):
    '''Get the mean of y.'''
    ymean = []
    for i in y:
        try:
            mm = ceil(stats.mean(i))
        except stats.StatisticsError as exc:
            ymean.append(0)
        else:
            ymean.append(mm)
    return ymean


def get_median(y):
    '''Get the median of y.'''
    ymedian = []
    for i in y:
        try:
            mm = ceil(stats.median(i))
        except stats.StatisticsError as exc:
            ymedian.append(0)
        else:
            ymedian.append(mm)
    return ymedian


def get_mode(y):
    '''Get the mode of y.'''
    ymode = []
    for i in y:
        try:
            mm = stats.mode(i)
        except stats.StatisticsError as exc:
            ymode.append(0)
        else:
            ymode.append(mm)
    return ymode


def get_stdev(y):
    '''Get the standard deviation of y.'''
    ystdev = []
    for i in y:
        if i == []:
            i = 0
        ystdev.append(np.std(i).item())
    return ystdev


def get_actual_mean_median_mode_stdev(pmus):
    ymus = list(pmus)
    mean = get_mean(ymus)
    median = get_median(ymus)
    mode = get_mode(ymus)
    stdev = get_stdev(ymus)
    return ymus, mean, median, mode, stdev


def config_compare_subplot(fig, ax, x):
    # Configure Y-axis
    ax.set_ylabel('Mean of the $\mu$ Estimates, $\mu_{mean}$',
                  fontweight='bold', fontsize=fontsize_axislabel)
    ax.tick_params(axis='both', which='both', labelsize=fontsize_axis)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(2))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(1))
    # Configure X-axis
    ax.set_xlabel('Dates', labelpad=10, fontweight='bold',
                  fontsize=fontsize_axislabel)
    fig.autofmt_xdate()  # rotate and align the tick labels so they look better. From https://matplotlib.org/3.1.3/gallery/recipes/common_date_problems.html
    ax.xaxis.set_major_locator(
        mdates.WeekdayLocator(byweekday=SU, interval=1, tz=None)
    )
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d %a", tz=None))
    ax.xaxis.set_minor_locator(
        mdates.WeekdayLocator(byweekday=(MO, TU, WE, TH, FR, SA),
                              interval=1, tz=None)
    )
    xmax = max(x)
    xmin = min(x)
    start = xmin - timedelta(days=2.)
    end = xmax + timedelta(days=2.)
    ax.set_xlim([start, end])
    # Draw Grid
    ax.grid(linestyle=':', color='grey')
    # Rotate xaxis minorticklabels - manual control
    plt.setp(ax.get_xminorticklabels(), rotation=90, ha='center')
    plt.setp(ax.get_xmajorticklabels(), rotation=90, ha='center')
    return fig, ax


def add_circuit_breaker(axes):
    cb_start = mdates.date2num(convert_str_to_datetime('07/04/2020'))
    cb_end = mdates.date2num(convert_str_to_datetime('01/06/2020'))
    width = cb_end - cb_start
    rect0 = Rectangle((cb_start, -1), width, 50, linewidth=0,
                      edgecolor=None, facecolor='grey', alpha=0.4)
    rect1 = Rectangle((cb_start, -20), width, 1500, linewidth=0,
                      edgecolor=None, facecolor='grey', alpha=0.4)
    axes[0].add_patch(rect0)
    axes[1].add_patch(rect1)
    fontcolor = 'black'
    msg = 'Circuit Breaker'
    msgx = mdates.date2num(convert_str_to_datetime('01/05/2020'))
    # msgy = 14
    # axes[0].text(msgx, msgy, msg, ha="left", va="bottom", color=fontcolor,
    #              fontsize=fontsize_axislabel, alpha=0.5)
    msgy = 400
    axes[1].text(msgx, msgy, msg, ha="left", va="bottom", color=fontcolor,
                 fontsize=fontsize_axislabel, alpha=0.8)
    return axes


def main(source1, source2, source3):
    # Read in predicted mus results to get x-axis and y-axis data
    # xdates is x-axis data
    predicted_mus1, xdates1 = read_predicted_mus(source1)
    predicted_mus2, xdates2 = read_predicted_mus(source2)
    predicted_mus3, xdates3 = read_predicted_mus(source3)
    # Get y-axis data
    # In situation where no mus are predicted, zero value is assigned to
    # facilitated the viewing of predicted mus
    ymus1, ymus1_mean, *_ = get_actual_mean_median_mode_stdev(
        predicted_mus1.values())
    ymus2, ymus2_mean, *_ = get_actual_mean_median_mode_stdev(
        predicted_mus2.values())
    ymus3, ymus3_mean, *_ = get_actual_mean_median_mode_stdev(
        predicted_mus3.values())
    # Get number of daily predictions of each attempt
    ymus1_len = [len(y) for y in ymus1]
    ymus2_len = [len(y) for y in ymus2]
    ymus3_len = [len(y) for y in ymus3]
    # Construct plot
    fig, ax = plt.subplots(2, 1, sharex=True, sharey=False)
    # Convert Python datetime to matplotlib date number. (forms the x-axis)
    mpdates = [mdates.date2num(date) for date in xdates1]
    # Add circuit-breaker zone
    ax = add_circuit_breaker(ax)
    # Plot line kwargs
    props1 = dict(color='r', alpha=0.5, lw=2, marker='d')
    props2 = dict(color='g', alpha=0.5, lw=2, marker='x')
    props3 = dict(color='b', alpha=0.5, lw=2, marker='o')
    props4 = dict(color='grey', alpha=0.5, lw=2, marker='d')
    props5 = dict(color='grey', alpha=0.5, lw=2, marker='x')
    props6 = dict(color='grey', alpha=0.5, lw=2, marker='o')
    ax[0].plot(mpdates, ymus1_mean, label='Seed1', **props1)
    ax[0].plot(mpdates, ymus2_mean, label='Seed2', **props2)
    ax[0].plot(mpdates, ymus3_mean, label='Seed3', **props3)
    ax[1].plot(mpdates, ymus1_len, label='Seed1', **props4)
    ax[1].plot(mpdates, ymus2_len, label='Seed2', **props5)
    ax[1].plot(mpdates, ymus3_len, label='Seed3', **props6)
    # Configure Subplots
    config_compare_subplot(fig, ax[0], xdates1)
    config_compare_subplot(fig, ax[1], xdates1)
    ax[1].get_yaxis().set_label_text(label=r'No. of $\mu$ Estimates',
                                     fontweight='bold')
    ax[0].set_ylim([-1, 18])
    ax[1].set_ylim([-20, 500])
    ax[1].yaxis.set_major_locator(ticker.MultipleLocator(100))
    ax[1].yaxis.set_minor_locator(ticker.MultipleLocator(10))
    # Draw Legend
    ax[0].legend(loc='upper left', fontsize=fontsize_axislabel, ncol=1)
    ax[1].legend(loc='lower left', fontsize=fontsize_axislabel, ncol=1)
    # Configure Subplot layout
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
    repo_dir = Path(__file__).resolve().parents[2]
    results_data = repo_dir / '2_Results' / '1_statistical_μ'
    fname1 = results_data / 'nd_predicted_mus_a1_s91023483923_pc8458.csv'
    fname2 = results_data / 'nd_predicted_mus_a2_s66576402094_pc8458.csv'
    fname3 = results_data / 'nd_predicted_mus_a3_s343090589475_pc8458.csv'
    main(fname1, fname2, fname3)
    plt.show()
