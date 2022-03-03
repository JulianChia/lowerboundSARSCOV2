from datetime import datetime, timedelta, timezone
from pathlib import Path
import csv
import random
import statistics as stats

from scipy.stats import norm, gaussian_kde
import numpy as np

from debug_functions import prv, pri, npri

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
from matplotlib.dates import MO, TU, WE, TH, FR, SA, SU

fontsize_axis = 16
fontsize_axislabel = 20


class CaseInfo(object):
    '''Class to store each COVID19 case info.

    Required Argument:

    args - It's elements are
       0      Case Number,
       1-3    Announced Confirmed Date, Announced Discharged Date,
              Announced Death Date,
       4-6    Actual Confirmed Date, Actual Discharged Date, Actual Death Date,
       7-9    Exposure, Imported from, Travel Info,
       10-12  Nationality, Age, Sex,
       13-14  Admission Venue, Admission Date,
       15-17  Quarantine Order Date, Quarantine Order Venue, Referred by MOH,
       18-19  Cluster, Relationship to other case,
       20-22  1st Symptom(s), 1st Symptom(s) Date, 1st Symptom(s) Venue,
       23-25  2nd Symptom(s), 2nd Symptom(s) Date, 2nd Symptom(s) Venue,
       26-27  1st Treatment Date, 1st Treatment Venue,
       28-29  2nd Treatment Date, 2nd Treatment Venue,
       30-31  3rd Treatment Date, 3rd Treatment Venue,
       32-33  4th Treatment Date, 4th Treatment Venue,
       34-35  5th Treatment Date, 5th Treatment Venue,
       36-37  6th Treatment Date, 6th Treatment Venue,
       38     Residence,
       39     Places Visited Prior to Hospital Admission
    '''

    def __init__(self, *args):
        # pri( args, args )
        self.number = int(args[0])
        self.announced_dates = {'Confirmed': self._datetimeit(args[1]),
                                'Discharged': self._datetimeit(args[2]),
                                'Deceased': self._datetimeit(args[3])
                                }
        self.actual_dates = {'Confirmed': self._datetimeit(args[4]),
                             'Discharged': self._datetimeit(args[5]),
                             'Deceased': self._datetimeit(args[6])
                             }
        self.exposures = {'Type': args[7], 'Imported From': args[8],
                          'Travel Info': args[9]}
        try:
            self.info = {'Nationality': args[10], 'Age': float(args[11]),
                         'Sex': args[12], 'Cluster': args[18],
                         'Links': args[19]
                         }
        except ValueError as exc:
            # print( f'Cases no. {self.number} {exc}' )
            # if 'could not convert string to float' in exc:
            self.info = {'Nationality': args[10], 'Age': None, 'Sex': args[12],
                         'Cluster': args[18], 'Links': args[19]
                         }
        self.symptoms = {1: (args[20], self._datetimeit(args[21]), args[22]),  # 1st : (symptoms, date, venue)
                         # 2nd : (symptoms, date, venue)
                         2: (args[23], self._datetimeit(args[24]), args[25])
                         }
        self.treatments = {1:  (self._datetimeit(args[26]), args[27]),  # 1st : (date, venue)
                           2:  (self._datetimeit(args[28]), args[29]),  # 2nd : (date, venue)
                           3:  (self._datetimeit(args[30]), args[31]),  # 3rd : (date, venue)
                           4:  (self._datetimeit(args[32]), args[33]),  # 4th : (date, venue)
                           5:  (self._datetimeit(args[34]), args[35]),  # 5th : (date, venue)
                           6:  (self._datetimeit(args[36]), args[37])   # 6th : (date, venue)
                           }
        self.quarantine_order = {'Venue': args[16], 'Date': args[15]}
        self.referred_by_MOH = args[17]
        self.admission = {'Venue': args[13],
                          'Date': self._datetimeit(args[14])}
        self.residence = args[38]
        self.places_visited = args[39]

    def _datetimeit(self, string_date):
        if string_date:
            try:
                date = datetime.strptime(string_date, "%d/%m/%Y")
            except IndexError:
                print(
                    f'\n    WARNING!!! Missing date in file "{source}" at line {line}. This data is excluded.\n')
            else:
                # print(date)
                datewtz = date.replace(tzinfo=timezone.utc)
            # print(f'\datewtz -> {datewtz}')
            return datewtz
        else:
            return string_date


def get_MOH_data(csvfile):
    '''Extract all data from COVID19 csv file. '''
    data = {}
    n = 1
    with open(csvfile, mode='r', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)  # Skip header row
        for row in reader:
            # print(f'row -> {type(row)} {row}')
            case = CaseInfo(*row)
            data[n] = case
            n += 1
    return data


def get_imported_covid19_arrival(data):

    def get_probable_days(mu, sigma):

        def quantile_function(mu, sigma):
            '''Returns the x of a Cumulative Distribution Function (CDF) for a known
            probability y=P(x). Here, 0 < P(x) <= 1 and is assigned randomly.

               mu    == mean, median or mode
               sigma == standard deviation

               A quantile function is also called the percent-point function
               or inverse cumulative distribution function.
            '''
            # y = random.random()
            y = random.uniform(0., 1.0)
            x = norm.ppf(y, loc=mu, scale=sigma)
            return x, y

        while True:
            days, probability = quantile_function(mu, sigma)
            if days > 0.:
                break
        return days, probability

    def _get_xy(time_axis, data):
        '''Return a dictionary storing each date in time_axis and the number of
        its occurrence in data.'''
        xy = {}
        for i in time_axis:
            cases = []
            for date, case in data:
                if i.date() == date.date():
                    cases.append(case)
            # print(i, len(cases))
            # print(i,  len(cases), cases)
            xy[i] = cases
        return xy

    def get_time_axis_of_event(all_dates, interval_days=1):
        dstart = min(all_dates)
        dend = max(all_dates) + timedelta(days=1)
        prv('dstart', dstart)
        prv('dend', dend)
        numdays = (dend - dstart).days
        prv('numdays', numdays)
        date_list = [dstart + timedelta(days=x) for x in range(numdays)]
        pri('date_list', date_list)
        return date_list

    def kde(x, x_grid, bandwidth=0.2, **kwargs):
        """Gaussian Kernel Density Estimation with Scipy"""
        kde = gaussian_kde(x, bw_method=bandwidth / x.std(ddof=1), **kwargs)
        return kde.evaluate(x_grid)

    def generate_rand_from_pdf(pdf, x_grid, size=1000, seed=datetime.now()):
        cdf = np.cumsum(pdf)
        cdf = cdf / cdf[-1]
        rng = np.random.default_rng(seed)
        values = rng.uniform(size=size)
        value_bins = np.searchsorted(cdf, values)
        random_from_cdf = x_grid[value_bins]
        return random_from_cdf

    # Get Imported COVID19 Case data from MOH data
    imported = [
        (v.number, v.announced_dates['Confirmed'], v.exposures['Travel Info'],
         v.admission['Date'])
        for v in data.values() if 'Import' in v.exposures['Type']
    ]
    pri('imported', imported)
    # Identify which COVID19 case has and not have arrival-date info.
    # - 1. Store them separate list with tuple elements of (date, case no.)
    # - 2. Also store their lag = confirm_date - arrival_date in
    #      another list.
    with_arrival_dates = []
    without_arrival_dates = []
    confirmation = []
    for i in imported:
        case = i[0]          # Case number
        confirm_date = i[1]  # Confirmed Date
        info = i[2]          # Arrival Info
        if info:
            arrival_date = convert_str_to_datetime(info[0:10])
            with_arrival_dates.append((arrival_date, case))
            lag = (confirm_date - arrival_date).days
            confirmation.append(lag)
        else:
            without_arrival_dates.append((confirm_date, case))
    # pri( 'confirmation', confirmation )
    pri('with_arrival_dates', with_arrival_dates)
    pri('without_arrival_dates', without_arrival_dates)
    pri('confirmation', confirmation)

    # Analyse mean and stdev of the time to confirm Imported Cases:
    # mean & stdev of Local COVID19 cases lag times
    confirmation_mean = stats.mean(confirmation)
    confirmation_stdev = stats.stdev(confirmation)
    confirmation_min = min(confirmation)
    confirmation_max = max(confirmation)
    prv('confirmation_mean', confirmation_mean)
    prv('confirmation_stdev', confirmation_stdev)
    prv('confirmation_min', confirmation_min)
    prv('confirmation_max', confirmation_max)

    fig, ax = plt.subplots(1, 1)
    # Empirical- CDF of with_arrival_dates
    hist_props = dict(density=True, cumulative=True, align='left',
                      histtype='bar', rwidth=0.8,
                      # histtype='stepfilled',
                      orientation='vertical', alpha=0.5)
    nbins = len(range(min(confirmation), max(confirmation)))
    prv('nbins', nbins)
    ax.hist(confirmation, bins=nbins,
            label=f'Empirical: {len(confirmation)} cases with published Arrival Dates',
            **hist_props)

    # Estimate Confirmation period for without_arrival_dates using
    # SciPy Kernel Density Estimation and conformation period of
    # with_arrival_dates.
    ihist, ibins = np.histogram(confirmation, bins=nbins)
    x_grid = np.linspace(confirmation_min, confirmation_max, len(confirmation))
    kdepdf = kde(np.array(confirmation), x_grid, bandwidth=0.1)
    random_from_kde = generate_rand_from_pdf(kdepdf, x_grid, seed=48945234)
    npri('random_from_kde', random_from_kde)
    label = f'Gaussian Kernel Density Estimation for {len(without_arrival_dates)} cases w/o arrival dates'
    ax.hist(random_from_kde, bins=nbins, label=label, **hist_props)

    # Normal Distribution Estimation
    label = r'$\mu=%4.3f, \sigma=%4.3f$' % (confirmation_mean, confirmation_stdev)
    # text = 'Gaussian Cumulative Distribution Function'
    fit_props = dict(color='green', linestyle='--', linewidth=0., marker='o',
                     markersize=5, alpha=0.8)
    n_data_points = 200
    x = np.linspace(norm.ppf(0.0001, confirmation_mean, confirmation_stdev),
                    norm.ppf(0.9999, confirmation_mean, confirmation_stdev),
                    n_data_points)
    y = norm.cdf(x, confirmation_mean, confirmation_stdev)
    ax.plot(x, y,
            label=f'Normal Distribution for {len(without_arrival_dates)} cases w/o arrival dates: {label}',
            **fit_props)

    ax.tick_params(axis='both', which='both', labelsize=fontsize_axis)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax.set_ylim([0, 1.1])
    axis_label_dict = dict(fontsize=fontsize_axislabel, fontweight='bold')
    ax.get_yaxis().set_label_text(
        label='Cumulative Probability', **axis_label_dict)
    ax.get_xaxis().set_label_text(
        label='Imported COVID-19 confirmation period since arrival to Singapore (Days)',
        **axis_label_dict)
    ax.grid(linestyle=':', color='grey')
    ax.legend(loc='upper left', fontsize=fontsize_axislabel)
    msg = 'Arrival to Singapore'
    msgx = -0.6
    msgy = 0.1
    ax.text(msgx, msgy, msg, ha="left", va="bottom", color='black',
            fontsize=fontsize_axislabel, alpha=0.8, rotation=90)
    ax.vlines(0, 0, 0.9, transform=ax.get_xaxis_transform(),
              colors='blue', linestyles='solid', alpha=0.8)
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

    plt.show()

    # Get real and probable Arrival time of Imported cases.
    imported_combined = [i for i in with_arrival_dates]
    # To estimate the arrival dates of those imported cases w/o arrival dates,
    # normal distribution is assumed.
    # arrival_date_estimated = [(cdate - timedelta(days=random_from_kde[n]), case) for n, (cdate, case) in enumerate(without_arrival_dates)]
    arrival_date_estimated = []
    for n, i in enumerate(without_arrival_dates):
        case = i[1]
        confirmeddate = i[0]  # Confirmed Date
        # days, probability = get_probable_days(mu, sigma)
        days = random_from_kde[n]
        arrival_date = confirmeddate - timedelta(days=days)
        arrival_date_estimated.append((arrival_date, case))
        imported_combined.append((arrival_date, case))
        # print( case, confirmeddate, days, probability, arrival_date )
    pri('##imported_combined', imported_combined)
    pri('##arrival_date_estimated', arrival_date_estimated)

    # Get time_axis of all Imported Cases
    arr_dates_real = {i for i, _ in with_arrival_dates}
    arr_dates_estimate = {i for i, _ in arrival_date_estimated}
    all_dates = arr_dates_real.union(arr_dates_estimate)
    time_axis_py = get_time_axis_of_event(all_dates)

    # Create dictionary {[time_axis] : [no. of cases]}
    with_arrival_dates = _get_xy(time_axis_py, with_arrival_dates)
    arrival_date_estimated = _get_xy(time_axis_py, arrival_date_estimated)
    imported_combined = _get_xy(time_axis_py, imported_combined)

    return with_arrival_dates, arrival_date_estimated, imported_combined


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
    # Plot Figure
    fig, ax = plt.subplots(3, 2, sharex=True, sharey=False)
    axes = (ax[0, 0], ax[1, 0], ax[2, 0], ax[0, 1], ax[1, 1], ax[2, 1])
    colors = ['C1', 'C6']
    ix = [mdates.date2num(date) for date in list(imported)]
    iy = [len(i) for i in imported.values()]
    pri('ix', ix)
    pri('iy', iy)

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
        axes[n].bar(ix, iy, label=r'Imported COVID-19 Arrival: Empirical & $KDE_{gauss}$',
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
        axes[n+3].bar(ix, iy, label=r'Imported COVID-19 Arrival: Empirical & $KDE_{gauss}$',
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
        ax.get_xaxis().set_label_text(label='Dates',
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
    cwd = Path(__file__).parent

    # 1. Get empirical imported COVID19 Cases and their arrival dates into
    #    Singapore.
    repo_dir = Path(__file__).resolve().parents[2]
    confirmed_cases = 'COVID19_Confirmed_Cases_Info_till_2020_05_27.csv'
    covid_data = repo_dir / '4_Empirical_Data' / confirmed_cases
    data = get_MOH_data(covid_data)
    with_arrival_dates, arrival_dates_estimated, imported_combined = \
        get_imported_covid19_arrival(data)
    pri('with_arrival_dates', with_arrival_dates)
    pri('arrival_dates_estimated', arrival_dates_estimated)
    pri('imported_combined', imported_combined)

    # 2. Get lower-bound Local SARS-CoV-2 infection trends
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
    print(CADS)
    print(WCADS)

    # Plot Results
    plot_SARSCoV2_curves(CADS, WCADS, imported_combined)
    # # seed = 91023483923
    # # seed = 66576402094
    # # seed = 343090589475
