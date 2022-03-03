import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
from matplotlib.dates import MO, TU, WE, TH, FR, SA, SU

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


def plot_SARSCoV2_curves(CADS, WCADS, imported):
    fontsize_axis = 16
    fontsize_axislabel = 20
    # Plot Figure
    fig, ax = plt.subplots(3, 2, sharex=True, sharey=False)
    axes = (ax[0, 0], ax[1, 0], ax[2, 0], ax[0, 1], ax[1, 1], ax[2, 1])
    colors = ['green', 'C6']
    ix = [mdates.date2num(date) for date in list(imported)]
    iy = list(imported.values())

    def add_zoom_region(ax, xdates, yp):
        # Inset axes of magnification
        # cb_start = mdates.date2num(convert_str_to_datetime('02/02/2020'))
        xo = mdates.date2num(convert_str_to_datetime('23/02/2020'))
        xe = mdates.date2num(convert_str_to_datetime('05/04/2020'))
        width = xe - xo
        bounds = [xo, 300, width, 1200]
        axins = ax.inset_axes(bounds, transform=ax.transData)
        # Plot COVID19 curves empirical and estimated
        axins.bar(xdates, yp, color=colors[1], alpha=0.5)
        axins.bar(ix, iy, color=colors[0], alpha=0.5)
        # Configure axins X-axis
        axins.xaxis.set_major_locator(
            mdates.WeekdayLocator(byweekday=SU, interval=1, tz=None))
        axins.xaxis.set_major_formatter(
            mdates.DateFormatter("%b %d %a", tz=None))
        axins.xaxis.set_minor_locator(
            mdates.WeekdayLocator(byweekday=(MO, TU, WE, TH, FR, SA),
                                  interval=1, tz=None))
        axins.set_xlim([xo, xe])
        # Configure ticks and border
        axins.tick_params(axis='both', which='both', labelsize=12,
                          color='blue', labelcolor='blue')
        for spine in axins.spines.values():
            spine.set_edgecolor('blue')
        # Configure Y-axis
        axins.yaxis.set_major_locator(ticker.MultipleLocator(5))
        axins.yaxis.set_minor_locator(ticker.MultipleLocator(1))
        axins.set_ylim([0, 50])
        # axins.yaxis.set_label_position("right")
        # axins.yaxis.tick_right()
        # Draw gridlines
        axins.grid(linestyle=':', color='blue')
        # Hide x-axis labels
        plt.setp(axins.get_xticklabels(), visible=False)
        # Draw boxs and connetion lines between ax and axins.
        ax.indicate_inset_zoom(axins, edgecolor="blue")

    # Plot subplots 1, 2, 3
    for n, (seed, v) in enumerate(CADS.items()):
        # Add circuit-breaker zone
        print()
        x_dates = [mdates.date2num(date) for date in v['pdates']]
        y_sarscov2 = v['inf_days_ncases']
        score = v['score'][1]
        # Plot data
        axes[n].bar(x_dates, y_sarscov2,
                    label=f'Local SARS-CoV-2: {seed}, CAD={score}',
                    color=colors[1], alpha=0.5)
        axes[n].bar(ix, iy, label='Imported COVID-19: Empirical',
                    color=colors[0], alpha=0.5)
        # Add Zoom Insert
        add_zoom_region(axes[n], x_dates, y_sarscov2)
    for n, (seed, v) in enumerate(WCADS.items()):
        # Add circuit-breaker zone
        print()
        x_dates = [mdates.date2num(date) for date in v['pdates']]
        y_sarscov2 = v['inf_days_ncases']
        score = round(float(v['score'][1]), 2)
        # Plot data
        axes[n+3].bar(x_dates, y_sarscov2,
                      label=f'Local SARS-CoV-2: {seed}, WCAD={score}',
                      color=colors[1], alpha=0.5)
        axes[n+3].bar(ix, iy, label='Imported COVID-19: Empirical',
                      color=colors[0], alpha=0.5)
        # Add Zoom Insert
        add_zoom_region(axes[n+3], x_dates, y_sarscov2)
        axes[n+3].yaxis.set_ticklabels([])
    for ax in axes:
        # Configure Y-axis
        # axes[n].get_yaxis().set_label_text(
        #     label='Number of Daily Cases',
        #     fontsize=fontsize_axislabel,
        #     fontweight='bold')
        ax.tick_params(axis='both', which='both', labelsize=fontsize_axis)
        ax.yaxis.set_major_locator(ticker.MultipleLocator(100))
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(20))
        ax.set_ylim([0, 2150])
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
        start = convert_str_to_datetime('16/02/2020')
        start = start - timedelta(days=3.)
        end = convert_str_to_datetime('19/04/2020')
        end = end + timedelta(days=3.)
        ax.set_xlim([start, end])
        # Draw Grid
        ax.grid(linestyle=':', color='grey')
        # Draw Legend
        ax.legend(loc='upper left', fontsize=fontsize_axislabel)

    # # Configure Subplot layput
    top = 0.980  # the top of the subplots of the figure
    bottom = 0.110  # the bottom of the subplots of the figure
    left = 0.060  # the left side of the subplots of the figure
    right = 0.985  # the right side of the subplots of the figure
    hspace = 0.050  # the amount of height reserved for space between subplots,
    # expressed as a fraction of the average axis height
    wspace = 0.050  # the amount of width reserved for space between subplots,
    # expressed as a fraction of the average axis width
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top,
                        wspace=wspace, hspace=hspace)
    # Rotate xaxis minorticklabels
    angle = 90
    plt.setp(axes[2].get_xminorticklabels(), rotation=angle, ha='center')  # manual control
    plt.setp(axes[2].get_xmajorticklabels(), rotation=angle, ha='center')  # manual control
    plt.setp(axes[5].get_xminorticklabels(), rotation=angle, ha='center')  # manual control
    plt.setp(axes[5].get_xmajorticklabels(), rotation=angle, ha='center')  # manual control
    axes[1].get_yaxis().set_label_text(label='Number of Daily Cases',
                                       fontsize=fontsize_axislabel,
                                       fontweight='bold')
    plt.show()


if __name__ == "__main__":
    # Source files with predicted mus
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
    # Get empirical imported COVID19 Cases
    covid19_csvfile = 'COVID19_epidemic_trends.csv'
    covid_data = repo_dir / '4_Empirical_Data' / covid19_csvfile
    imported, _ = get_daily_cases(covid_data)
    # Plot Results
    plot_SARSCoV2_curves(CADS, WCADS, imported)
    # seed = 91023483923
    # seed = 66576402094
    # seed = 343090589475
