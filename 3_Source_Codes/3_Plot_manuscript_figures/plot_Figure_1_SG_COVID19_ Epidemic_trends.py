#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Program to plot Singapore's imported and local COVID-19 epidemic curves.
'''
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
from matplotlib.dates import MO, TU, WE, TH, FR, SA, SU
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from matplotlib.patches import Rectangle

from itertools import count
from datetime import datetime, timedelta, timezone
import csv
from pathlib import Path

fontsize_axis = 16
fontsize_axislabel = 20


def _get_daily_cases(csvfile):
    '''Extract the number of Daily Imported and Local Confirmed COVID19 cases
       from "csvfile" and return corresponding dictionaries.
    '''
    imported = {}
    local = {}
    with open(csvfile, mode='r', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)  # Skip header row
        for row in reader:
            if row:
                # print( f'row -> {type(row)} {row}')
                date = convert_str_to_datetime(row[0], dateformat='%Y-%m-%d')
                imported[date] = int(row[1])
                local[date] = int(row[2])
                # print( date, imported, local )

    index = count()
    for k, v in imported.items():
        print(f'{next(index)} {k}  {v}  {local[k]}')

    return imported, local


def convert_str_to_datetime(date, dateformat="%d/%m/%Y", tzinfo=timezone.utc):
    refdate = datetime.strptime(date, dateformat)
    refdate = refdate.replace(tzinfo=tzinfo)
    return refdate


def plot_covid19_epidemic_curves(imported, local):
    # Plot Figure
    fig, ax = plt.subplots(2, 1, sharex=True, sharey=False)
    # fig.suptitle('Singapore COVID-19 Epidemic Curves', fontsize=16,
    #              fontweight='bold')
    # Add circuit-breaker zone
    ax = add_circuit_breaker(ax)
    # Plot Empirical COVID19 data
    ix = [mdates.date2num(date) for date in list(imported.keys())]
    iy = list(imported.values())
    lx = [mdates.date2num(date) for date in list(local.keys())]
    ly = list(local.values())
    ilabel = 'Imported COVID19 cases'
    llabel = 'Local COVID19 cases'
    ax[0].bar(ix, iy, label=ilabel, color='green', alpha=0.5)
    ax[1].bar(lx, ly, label=llabel, color='C0', alpha=0.5)
    # Configure Y-axis
    ax[0].get_yaxis().set_label_text(label='Number of Daily Cases',
                                     fontsize=fontsize_axislabel,
                                     fontweight='bold')
    ax[0].tick_params(axis='both', which='both', labelsize=fontsize_axis)
    ax[0].yaxis.set_major_locator(ticker.MultipleLocator(5))
    ax[0].yaxis.set_minor_locator(ticker.MultipleLocator(1))
    ax[0].set_ylim([0, 50])

    ax[1].get_yaxis().set_label_text(label='Number of Daily Cases',
                                     fontsize=fontsize_axislabel,
                                     fontweight='bold')
    ax[1].tick_params(axis='both', which='both', labelsize=fontsize_axis)
    ax[1].yaxis.set_major_locator(ticker.MultipleLocator(100))
    ax[1].yaxis.set_minor_locator(ticker.MultipleLocator(20))
    ax[1].set_ylim([0, 1500])
    # Configure X-axis
    ax[1].get_xaxis().set_label_text(label='Press Release Dates',
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
    # Draw Grid
    ax[0].grid(linestyle=':', color='grey')
    ax[1].grid(linestyle=':', color='grey')
    # Draw Legend
    ax[0].legend(loc='upper right', fontsize=fontsize_axislabel)
    ax[1].legend(loc='upper right', fontsize=fontsize_axislabel)
    # Rotate xaxis minorticklabels
    angle = 90
    plt.setp(ax[1].get_xminorticklabels(), rotation=angle, ha='center')  # manual control
    plt.setp(ax[1].get_xmajorticklabels(), rotation=angle, ha='center')  # manual control
    # Configure Subplot layput
    top = 0.980  # the top of the subplots of the figure
    bottom = 0.130  # the bottom of the subplots of the figure
    left = 0.070  # the left side of the subplots of the figure
    right = 0.985  # the right side of the subplots of the figure
    hspace = 0.060  # the amount of height reserved for space between subplots,
    # expressed as a fraction of the average axis height
    wspace = 0.220  # the amount of width reserved for space between subplots,
    # expressed as a fraction of the average axis width
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top,
                        wspace=wspace, hspace=hspace)
    # Insert Zoomed in region
    axins = add_zoom_region(ax[1], lx, ly)
    return fig, ax


def add_circuit_breaker(axes):
    cb_start = mdates.date2num(convert_str_to_datetime('07/04/2020'))
    cb_end = mdates.date2num(convert_str_to_datetime('01/06/2020'))
    width = cb_end - cb_start

    rect0 = Rectangle((cb_start, 0), width, 50, linewidth=0, edgecolor=None,
                      facecolor='grey', alpha=0.4)
    rect1 = Rectangle((cb_start, 0), width, 1500, linewidth=0, edgecolor=None,
                      facecolor='grey', alpha=0.4)
    axes[0].add_patch(rect0)
    axes[1].add_patch(rect1)

    msg = 'Circuit Breaker'
    msgx = mdates.date2num(convert_str_to_datetime('26/04/2020'))
    msgy = 15
    fontcolor = 'black'
    axes[0].text(msgx, msgy, msg, ha="left", va="bottom", color=fontcolor,
                 fontsize=fontsize_axislabel, alpha=0.8)
    # msgy = 1300
    # axes[1].text(msgx, msgy, msg, ha="left", va="bottom", color=fontcolor,
    #              fontsize=16, alpha=0.5)
    return axes


def add_zoom_region(ax, lx, ly):
    # Inset axes of magnification
    xo = mdates.date2num(convert_str_to_datetime('02/02/2020'))
    xe = mdates.date2num(convert_str_to_datetime('29/03/2020'))
    width = xe - xo
    print('width = ', width)
    bounds = [xo, 300, width, 700]
    axins = ax.inset_axes(bounds, transform=ax.transData)
    # Plot dats
    axins.bar(lx[0:67], ly[0:67], color='C0', alpha=0.4)
    # Draw boxs and connection lines between ax and axins.
    mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="black")
    # Configure axins X-axis
    axins.xaxis.set_major_locator(
        mdates.WeekdayLocator(byweekday=SU, interval=1, tz=None))
    axins.xaxis.set_major_formatter(mdates.DateFormatter("%b %d %a", tz=None))
    axins.xaxis.set_minor_locator(
        mdates.WeekdayLocator(byweekday=(MO, TU, WE, TH, FR, SA),
                              interval=1, tz=None))
    x1 = mdates.date2num(convert_str_to_datetime('02/02/2020'))
    x2 = mdates.date2num(convert_str_to_datetime('29/03/2020'))
    axins.set_xlim([x1, x2])
    # Configure Y-axis
    axins.tick_params(axis='both', which='both', labelsize=14)
    axins.yaxis.set_major_locator(ticker.MultipleLocator(5))
    axins.yaxis.set_minor_locator(ticker.MultipleLocator(1))
    axins.set_ylim([0, 35])
    # Draw gridlines
    axins.grid(linestyle=':', color='blue')
    # Hide x-axis labels
    plt.setp(axins.get_xticklabels(), visible=False)
    return axins


def main():
    repo_dir = Path(__file__).resolve().parents[2]
    covid19_csvfile = 'COVID19_epidemic_trends.csv'
    covid_data = repo_dir / '4_Empirical_Data' / covid19_csvfile
    imported, local = _get_daily_cases(covid_data)
    fig, ax = plot_covid19_epidemic_curves(imported, local)


if __name__ == "__main__":
    main()
    plt.show()
