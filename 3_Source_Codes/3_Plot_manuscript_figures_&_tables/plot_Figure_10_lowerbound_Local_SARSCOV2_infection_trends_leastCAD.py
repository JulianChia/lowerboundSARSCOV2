import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
from matplotlib.dates import MO, TU, WE, TH, FR, SA, SU
from matplotlib.patches import Rectangle

from pathlib import Path
from datetime import datetime, timedelta, timezone
import csv


def get_daily_cases(csvfile):
    '''Extract the number of Daily Imported and Local Confirmed COVID19
    cases from "csvfile" and return corresponding dictionaries.
    '''
    imported = {}
    local = {}
    with open(csvfile, mode='r', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)  # Skip header row
        for row in reader:
            if row:
                # print( f'row -> {type(row)} {row}')
                date = convert_str_to_datetime(row[0],
                                               dateformat='%Y-%m-%d')
                imported[date] = int(row[1])
                local[date] = int(row[2])
                # print( date, imported, local )
    # for k, v in imported.items():
    #    print(f'{next(index)} {k}  {v}  {local[k]}')
    return imported, local


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


def plot_SARSCoV2_curves(CADS):
    fontsize_axis = 16
    fontsize_axislabel = 20

    def add_circuit_breaker(ax, subplot):
        cb_ann_x = mdates.date2num(convert_str_to_datetime('03/04/2020'))
        ax.vlines(cb_ann_x, 0, 1.0, transform=ax.get_xaxis_transform(),
                  colors='orange', linestyles='solid', alpha=0.8)
        cb_start = mdates.date2num(convert_str_to_datetime('07/04/2020'))
        cb_end = mdates.date2num(convert_str_to_datetime('01/06/2020'))
        width = cb_end - cb_start
        # Boundary
        rect = Rectangle((cb_start, 0), width, 1500, linewidth=0,
                         edgecolor=None, facecolor='grey', alpha=0.4)
        ax.add_patch(rect)
        if subplot == 1:
            # Text
            msg = 'Circuit Breaker'
            fontcolor = 'black'
            msgx = mdates.date2num(convert_str_to_datetime('03/05/2020'))
            msgy = 1100
            ax.text(msgx, msgy, msg, ha="left", va="bottom", color=fontcolor,
                    fontsize=fontsize_axislabel, alpha=0.8)
            msgx = mdates.date2num(convert_str_to_datetime('01/04/2020'))
            msgy = 400
            msg = r'$3^{rd}$ Apr: Circuit Breaker Announced.'
            ax.text(msgx, msgy, msg, ha="left", va="bottom", color="blue",
                    fontsize=12, rotation=90, alpha=0.8)

    def add_zoom_region(ax, xdates, ye, yp, colors):
        # Inset axes of magnification
        # cb_start = mdates.date2num(convert_str_to_datetime('02/02/2020'))
        xo = mdates.date2num(convert_str_to_datetime('19/01/2020'))
        xe = mdates.date2num(convert_str_to_datetime('29/03/2020'))
        width = xe - xo
        bounds = [xo, 300, width, 800]
        axins = ax.inset_axes(bounds, transform=ax.transData)
        # Plot COVID19 curves empirical and estimated
        axins.bar(xdates, ye, color=colors[0], alpha=0.5)
        axins.bar(xdates, yp, color=colors[1], alpha=0.5)

        # Configure axins X-axis
        axins.xaxis.set_major_locator(
            mdates.WeekdayLocator(byweekday=SU, interval=1, tz=None))
        axins.xaxis.set_major_formatter(
            mdates.DateFormatter("%b %d %a", tz=None))
        axins.xaxis.set_minor_locator(
            mdates.WeekdayLocator(byweekday=(MO, TU, WE, TH, FR, SA),
                                  interval=1, tz=None))
        axins.set_xlim([xo, xe])
        # Configure Y-axis
        axins.tick_params(axis='both', which='both', labelsize=12)
        axins.yaxis.set_major_locator(ticker.MultipleLocator(5))
        axins.yaxis.set_minor_locator(ticker.MultipleLocator(1))
        axins.set_ylim([0, 40])
        # axins.yaxis.set_label_position("right")
        # axins.yaxis.tick_right()
        # Draw gridlines
        axins.grid(linestyle=':', color='blue')
        # Hide x-axis labels
        plt.setp(axins.get_xticklabels(), visible=False)
        # Draw boxs and connetion lines between ax and axins.
        ax.indicate_inset_zoom(axins, edgecolor="black")

    # Plot Figure
    fig, ax = plt.subplots(3, 1, sharex=True, sharey=False)
    axes = (ax[0], ax[1], ax[2])
    colors = ['C0', 'C6']
    xmin = []
    xmax = []
    # Plot subplots 1, 2, 3
    for n, (seed, v) in enumerate(CADS.items()):
        # Add circuit-breaker zone
        add_circuit_breaker(ax[n], n)
        print()
        x_dates = [mdates.date2num(date) for date in v['pdates']]
        xmin.append(v['pdates'][0])
        xmax.append(v['pdates'][-1])
        y_emp_covid19 = v['plncases']
        y_sarscov2 = v['inf_days_ncases']
        score = v['score'][1]
        # Plot data
        axes[n].bar(x_dates, y_emp_covid19, label='Local COVID-19: Empirical',
                    color=colors[0], alpha=0.5)
        axes[n].bar(x_dates, y_sarscov2,
                    label=f'Local SARS-CoV-2: {seed}, CAD={score}',
                    color=colors[1], alpha=0.5)
        # Configure Y-axis
        axes[n].get_yaxis().set_label_text(
            label='Number of Daily Cases',
            fontsize=fontsize_axislabel,
            fontweight='bold')
        axes[n].tick_params(axis='both', which='both', labelsize=fontsize_axis)
        axes[n].yaxis.set_major_locator(ticker.MultipleLocator(100))
        axes[n].yaxis.set_minor_locator(ticker.MultipleLocator(20))
        axes[n].set_ylim([0, 1500])
        # Draw Grid
        ax[n].grid(linestyle=':', color='grey')
        # Draw Legend
        ax[n].legend(loc='upper right', fontsize=fontsize_axislabel)
        # Add Zoom Insert
        add_zoom_region(ax[n], x_dates, y_emp_covid19, y_sarscov2, colors)

    # Configure X-axis
    ax[2].get_xaxis().set_label_text(label='Press Release Dates',
                                     fontsize=fontsize_axislabel,
                                     fontweight='bold')
    fig.autofmt_xdate()  # rotate and align the tick labels so they look better. From https://matplotlib.org/3.1.3/gallery/recipes/common_date_problems.html
    ax[2].xaxis.set_major_locator(
        mdates.WeekdayLocator(byweekday=SU, interval=1, tz=None))
    ax[2].xaxis.set_major_formatter(mdates.DateFormatter("%b %d %a", tz=None))
    ax[2].xaxis.set_minor_locator(
        mdates.WeekdayLocator(byweekday=(MO, TU, WE, TH, FR, SA),
                              interval=1, tz=None))
    #ax.xaxis.set_minor_formatter(mdates.DateFormatter("%a", tz=None))
    # ax.set_xlim(auto=True)
    # start = convert_str_to_datetime('23/01/2020')
    start = convert_str_to_datetime('19/01/2020')
    start = start - timedelta(days=5.)
    # end = convert_str_to_datetime('21/08/2020')
    end = convert_str_to_datetime('18/08/2020')
    end = end + timedelta(days=2.)
    ax[2].set_xlim([start, end])
    # Rotate xaxis minorticklabels
    angle = 90
    plt.setp(ax[2].get_xminorticklabels(), rotation=angle, ha='center')  # manual control
    plt.setp(ax[2].get_xmajorticklabels(), rotation=angle, ha='center')  # manual control
    ax[0].get_yaxis().set_label_text(label=' ', fontsize=fontsize_axislabel,
                                     fontweight='bold')
    ax[2].get_yaxis().set_label_text(label=' ', fontsize=fontsize_axislabel,
                                     fontweight='bold')
    # Configure Subplot layput
    top = 0.980  # the top of the subplots of the figure
    bottom = 0.110  # the bottom of the subplots of the figure
    left = 0.060  # the left side of the subplots of the figure
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
    seeds = ('Seed1', 'Seed2', 'Seed3')
    CADS = {seeds[n]: read_CAD_WCAD(source=i) for n, i in enumerate(cadfiles)}
    # Plot Results
    plot_SARSCoV2_curves(CADS)
    # seed = 91023483923
    # seed = 66576402094
    # seed = 343090589475
