#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Program to plot the COVID-19 Confirmation Period (CCP) of Singapore's local
COVID-19 cases, their cumulative distribution, and the range of constant
average daily CCP, i.e. mu_c, for normal and truncated normal distributions
models.
'''
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as patches
import numpy as np
import statistics
from datetime import datetime, timezone
from pathlib import Path
from scipy.stats import norm
import csv


class CaseInfo(object):
    '''Class to store each COVID19 case info.

    Required Argument:

    args - It's elements are
       0      Case Number,
       1-3    Announced Confirmed Date, Announced Discharged Date, Announced
              Death Date,
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
        self.number = int(args[0])
        self.announced_dates = {'Confirmed': self._datetimeit(args[1]),
                                'Discharged': self._datetimeit(args[2]),
                                'Deceased': self._datetimeit(args[3])
                                }
        self.actual_dates = {'Confirmed': self._datetimeit(args[4]),
                             'Discharged': self._datetimeit(args[5]),
                             'Deceased': self._datetimeit(args[6])
                             }
        self.exposures = {'Type': args[7],
                          'Imported From': args[8],
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
        self.symptoms = {1: (args[20], self._datetimeit(args[21]), args[22]),  # 1st:(symptoms, date, venue)
                         # 2nd:(symptoms, date, venue)
                         2: (args[23], self._datetimeit(args[24]), args[25])
                         }
        self.treatments = {1: (self._datetimeit(args[26]), args[27]),  # 1st:(date, venue)
                           2: (self._datetimeit(args[28]), args[29]),  # 2nd:(date, venue)
                           3: (self._datetimeit(args[30]), args[31]),  # 3rd:(date, venue)
                           4: (self._datetimeit(args[32]), args[33]),  # 4th:(date, venue)
                           5: (self._datetimeit(args[34]), args[35]),  # 5th:(date, venue)
                           6: (self._datetimeit(args[36]), args[37])  # 6th:(date, venue)
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
                datewtz = date.replace(tzinfo=timezone.utc)
            return datewtz
        else:
            return string_date


def get_MOH_COVID19_data(csvfile):
    '''Extract all data from COVID19 csv file. '''
    data = {}
    n = 1
    with open(csvfile, mode='r', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)  # Skip header row
        for row in reader:
            # print( f'row -> {type(row)} {row}')
            case = CaseInfo(*row)
            data[n] = case
            n += 1
    return data


def prv(name, arg):
    print(f'\n{name} -> {type(arg)} -> {arg}')


def pri(name, arg):
    print(f'\n{name} -> {len(arg)} -> {arg}')


def convert_str_to_pydate(date, dateformat="%d/%m/%Y", tzinfo=timezone.utc):
    refdate = datetime.strptime(date, dateformat)
    refdate = refdate.replace(tzinfo=tzinfo)
    return refdate


def _get_dayncases(days_axis, data):
    xy = {}
    for i in days_axis:
        cases = []
        for days, case in data:
            if i == days:
                cases.append(case)
        # print( i, cases)
        xy[i] = cases
    return xy


def get_known_covid19_Confirmation_Period(data):
    events = [convert_str_to_pydate('19/01/2020'),  # Life Church Missions & Missions Singapore service
              # Tour guide first meet Wuhan Tour group and infect husband
              convert_str_to_pydate('22/01/2020'),
              # Yong Thai Hang staff meets Wuhan Tour group & infect husband
              convert_str_to_pydate('23/01/2020'),
              convert_str_to_pydate('25/01/2020'),  # CNY family gathering at Mei Hwan Drive
              convert_str_to_pydate('29/01/2020'),  # Grace Assembly of God staff devotion meeting
              # Grace Assembly of God staff devotion meeting & prayer meeting with attendees
              convert_str_to_pydate('05/02/2020'),
              convert_str_to_pydate('15/02/2020')   # Jurong Safra
              ]
    clusters_cases = ([31,  33,  38,  83,  91],
                      [24,  25],
                      [19,  20,  21,  27],
                      [66,  68,  70,  71,  80,  84,  88],
                      [48,  49,  51],
                      [53,  54,  57,  60,  61,  62,  63,  67,  73],
                      [107, 112, 114, 115, 116, 117, 118, 120, 121,
                       122, 125, 127, 128, 129, 130, 131, 133, 134, 137]
                      )
    pri('clusters_cases', clusters_cases)

    confirmed = []
    for i, cluster in enumerate(clusters_cases):
        c = get_refcluster_data(data, cluster, events[i])
        confirmed.append(c)
    pri('confirmed', confirmed)

    cdays = [j[0].days for i in confirmed for j in i]
    # pri('cdays', cdays)

    first = min(cdays)
    if first >= 0:
        first = 0
    last = max(cdays) + 1
    daysline = [i for i in range(first, last, 1)]
    # pri('daysline', daysline)

    clist = [(j[0].days, j[1]) for i in confirmed for j in i]
    # pri('clist', clist)

    # Store results in dictionaries {ay number:[cases]}
    confirmed_dict = _get_dayncases(daysline, clist)
    pri('confirmed_dict', confirmed_dict)

    return confirmed_dict, clusters_cases, confirmed


def get_refcluster_data(data, cases, eventdate):
    confirmed = [(data[i].announced_dates['Confirmed']-eventdate, data[i].number)
                 for i in cases if data[i].announced_dates['Confirmed']]
    # pri('confirmed', confirmed )
    return confirmed


def plot_CCPs_and_constant_average_daily_CCPs(confirmed, clusters_cases, con):
    fontsize_axis = 16
    fontsize_axislabel = 20

    def add_cluster_labels(ax, data_type, clusters, bbox_props, fontsize=8):
        x = list(data_type.keys())
        y = list(data_type.values())
        pri('x', x)
        pri('y', y)
        ann_list = []
        fontsize = fontsize

        def choose_color(case, clusters, bbox_props):
            for i, cluster in enumerate(clusters):
                if case in cluster:
                    return bbox_props[i]

        # Annotate Stemline
        for i, cases in enumerate(y):
            count = 1
            for case in cases:
                bbox_properties = choose_color(case, clusters, bbox_props)
                ann = ax.annotate(
                    f'{case}', xy=(x[i], count), xycoords='data', picker=True,
                    color='blue', rotation=0, ha='center', va='center',
                    bbox=bbox_properties, fontsize=fontsize)
                ann_list.append(ann)
                count += 1
        return ann_list

    # 1. Create figure and subplots
    fig = plt.figure()
    ax1 = plt.subplot2grid((3, 1), (0, 0))
    ax2 = plt.subplot2grid((3, 1), (1, 0), sharex=ax1)
    ax3 = plt.subplot2grid((3, 1), (2, 0), sharex=ax1)

    # 2. Draw stemlines in subplots 1
    x = list(confirmed.keys())
    cy = [len(i) for i in confirmed.values()]
    ccolor = '#0000ff'
    csp = ax1.stem(x, cy, markerfmt=None, linefmt="-", basefmt="k-",
                   use_line_collection=True)
    csp[0].set_marker(None)
    csp[1].set_color(ccolor)
    text = '- Documented COVID-19\n  case number in blue.'
    ax1.text(31, 2.5, text, ha="left", va="baseline", color="blue",
             fontsize=fontsize_axislabel-2, alpha=0.6,)
    stemplots = []
    stemplots.append(csp)
    pri('stemplots', stemplots)
    # 2.1 Add Annotations into Stemplots
    pad = 0.2
    c1 = '#a6ff85'   # light green
    c2 = 'cyan'      # cyanish
    c3 = '#b4aaff'   # light purple
    c4 = '#ff85a6'   # light red
    c5 = '#ffb285'   # light orange
    c6 = '#85c6ff'   # light blue
    c7 = '#fffb8b'   # beige yellow
    ccs = (c1, c2, c3, c4, c5, c6, c7)
    width = 4
    height = 1
    bbox_props = []
    custom_rect = []
    for cc in ccs:
        bbox_props.append(
            dict(boxstyle=f'round,pad={pad}', fc=cc, alpha=1., ec=cc)
        )
        custom_rect.append(
            patches.Rectangle((0, 0), width, height, fc=cc, alpha=1., ec=cc)
        )
    add_cluster_labels(ax1, confirmed, clusters_cases, bbox_props, fontsize=13)
    # 2.2 Change Legend labels and handles
    custom_label = ['Infected on Jan 19',
                    'Infected on Jan 22',
                    'Infected on Jan 23',
                    'Infected on Jan 25',
                    'Infected on Jan 29',
                    'Infected on Feb 5 ',
                    'Infected on Feb 15 ']
    ax1.legend(custom_rect, custom_label, loc='upper right',
               fontsize=fontsize_axislabel, ncol=1,)

    # 3. Draw Empirical Cumulative Probability Distribution
    # subplot 2 - emprical
    con_days = [j[0].days for i in con for j in i]
    con_bins = [i for i in range(0, max(con_days)+1, 1)]
    hist_props = dict(density=True, cumulative=True, align='mid',
                      # histtype='bar', rwidth=0.8,
                      histtype='stepfilled', orientation='vertical', alpha=0.2)
    con_hist = ax2.hist(con_days, bins=con_bins, color=ccolor,
                        label='Empirical', **hist_props)
    days_mu = statistics.mean(con_days)
    days_sigma = statistics.pstdev(con_days)
    label = r'$\mu=%4.3f, \sigma=%4.3f$' % (days_mu, days_sigma)
    text = 'Gaussian Cumulative Distribution Function'
    n_data_points = 200
    fit_props = dict(color='black', linestyle='--', linewidth=0., marker='o',
                     markersize=3)
    # subplot 2 - normal distribution
    x = np.linspace(norm.ppf(0.0001, days_mu, days_sigma),
                    norm.ppf(0.9999, days_mu, days_sigma),
                    n_data_points)
    y = norm.cdf(x, days_mu, days_sigma)
    ax2.plot(x, y, label=label, **fit_props)
    ax2.text(17.3, 0.5, text, ha="right", va="baseline", color="black",
             fontsize=fontsize_axislabel, alpha=1.,)
    # Configure Legend
    ax2.get_legend_handles_labels()
    ax2.legend(prop={'size': fontsize_axislabel}, loc='upper left')

    # 4. Draw Modelled Normal Cumulative Probability Distribution in subplot 3
    def _get_mu_sigma(mu_min_max, sigma_constants):
        '''Function to define the average number of days taken to confirm
        local covid19 cases and their corresponding standard deviations.

        Arguments:
         mu_min_max - min. and max. average number of days taken to locally
                      confirm COVID19 after SARS-CoV2 infection.
         sigma_constants - curve fitting constants (a,b,c) for a quadratic
                           mu-sigma relationship.

        Return:
         mu_sigma - a dictionary of the mu-sigma relationship
                    key = a range of the average number of days taken to
                          locally confirm COVID19 after SARS-CoV2 infection
                    value = standard deviations of the key
        '''
        mu_min, mu_max = mu_min_max
        a, b, c = sigma_constants
        mu = [i for i in range(mu_min, mu_max+1, 1)]
        mu_sigma = {x: a*x**2+b*x+c for x in mu}
        pri('mu_sigma', mu_sigma)
        return mu_sigma

    mu_min_max = (1, 18)
    sigma_constants = (-0.008665, 0.483888, 0.0)
    mu_sigma = _get_mu_sigma(mu_min_max, sigma_constants)
    mu_sigma_probability = []
    left = 0
    right = np.inf
    for k, v in mu_sigma.items():
        mu = k
        sigma = v
        mu_sigma_label = r'$\mu_c=%2.1f, \sigma_c=%4.3f$' % (mu, sigma)
        # subplot 3
        x = np.linspace(norm.ppf(0.0001, mu, sigma),
                        norm.ppf(0.9999, mu, sigma),
                        n_data_points)
        y = norm.cdf(x, mu, sigma)
        mu_sigma_probability.append(ax3.plot(x, y, label=mu_sigma_label,
                                             linestyle='--', linewidth=0.,
                                             marker='o', markersize=3))
    # Configure Legend
    ax3.get_legend_handles_labels()
    ax3.legend(prop={'size': 14}, loc='center right', ncol=2)
    # Add text to show when SARS-CoV-2 infection started
    textmu_props = dict(ha="center", va="baseline", color="black",
                        fontsize=fontsize_axislabel-2, alpha=0.6,)
    textlegend_title_props = dict(ha="left", va="baseline", color="black",
                                  fontsize=fontsize_axislabel, alpha=1.0,)
    textlegend_props = dict(ha="left", va="baseline", color="black",
                            fontsize=fontsize_axis, alpha=0.6,)
    textmu1 = r'$\mu_c=1$'
    textmu18 = r'$\mu_c=18$'
    textlegend = 'Gaussian Cumulative Distribution'
    texteqn1 = r'$\sigma_c = a\mu_c^2 + b\mu_c + c$'
    texteqn2 = r'$a=%f, b=%f, c=%2.1f$' % (sigma_constants)
    ax3.text(0.75, 0.48, textmu1, rotation=80, **textmu_props)
    ax3.text(19.2, 0.48, textmu18, rotation=30, **textmu_props)
    ax3.text(27.5, 0.85, textlegend,  **textlegend_title_props)
    ax3.text(29.0, 0.14, texteqn1, **textlegend_props)
    ax3.text(29.0, 0.08, texteqn2, **textlegend_props)

    # 5. Configure subplots y-axis
    ax1.get_yaxis().set_label_text(label='No. of Cases', fontweight='bold',
                                   fontsize=fontsize_axislabel,)
    ax1.yaxis.set_major_locator(ticker.MultipleLocator(5))
    ax1.yaxis.set_minor_locator(ticker.MultipleLocator(1))
    ax1.set_ylim([-0.5, max(cy)+0.5])

    p_ax_list = [ax2, ax3]
    for ax in p_ax_list:
        ax.get_yaxis().set_label_text(label='Cumulative Probability',
                                      fontweight='bold',
                                      fontsize=fontsize_axislabel,)
        ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
        # ax.set_ylim(auto = True)
        ax.set_ylim([-0.05, 1.05])
        ax.tick_params(axis='both', which='both', labelsize=fontsize_axis)

    # 5. Configure plot x-axis
    ax1.xaxis.set_major_locator(ticker.MultipleLocator(7))
    ax1.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    # ax1.set_xlim( auto=True )
    ax1.set_xlim([-4, 40])

    for ax in p_ax_list:
        ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
        # ax.set_xlim(auto=True)
        ax.set_xlim([-4, 40])
    xlabel = 'Number of Days After SARS-CoV-2 Infection'
    ax3.get_xaxis().set_label_text(label=xlabel, fontweight='bold',
                                   fontsize=fontsize_axislabel,)

    # 6. Configure Subplots Grid
    ax_list = [ax1, ax2, ax3]
    for ax in ax_list:
        ax.grid(linestyle='--', color='grey')
        ax.grid(linestyle=':',  color='grey', which='minor')

    # 7. Configure Subplot layput
    top = 0.980  # the top of the subplots of the figure
    bottom = 0.050  # the bottom of the subplots of the figure
    left = 0.050  # the left side of the subplots of the figure
    right = 0.990  # the right side of the subplots of the figure
    hspace = 0.000  # the amount of height reserved for space between subplots,
    # expressed as a fraction of the average axis height
    wspace = 0.170  # the amount of width reserved for space between subplots,
    # expressed as a fraction of the average axis width
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top,
                        wspace=wspace, hspace=hspace)

    # 9. Add text to show when SARS-CoV-2 infection started
    ax_list = [ax2, ax3]
    stdict = dict(ha="center", va="baseline", color="black",
                  fontsize=fontsize_axislabel, rotation=90, alpha=0.6,)
    for i in ax_list:
        i.text(0.0, 0.05, 'SARS-CoV-2 Infection', **stdict)
    ax1.text(0.0, 0.5, 'SARS-CoV-2 Infection', **stdict)

    # 13. Hide subplot1,2,3 xmajorticklabels
    plt.setp(ax1.get_xmajorticklabels(), visible=False)  # manual control
    plt.setp(ax2.get_xmajorticklabels(), visible=False)  # manual control
    # plt.setp(ax3.get_xmajorticklabels(), visible=False)  # manual control

    # change majorticklabels fontsize
    ax1.tick_params(axis='both', labelsize=fontsize_axis)
    ax2.tick_params(axis='both', labelsize=fontsize_axis)
    # ax3.tick_params(axis='both', labelsize=fontsize_axis)

    return fig


def main():
    repo_dir = Path(__file__).resolve().parents[2]
    confirmed_cases = 'COVID19_Confirmed_Cases_Info_till_2020_05_27.csv'
    covid_data = repo_dir / '4_Empirical_Data' / confirmed_cases
    moh = get_MOH_COVID19_data(covid_data)
    ccp_data = get_known_covid19_Confirmation_Period(moh)
    plot_CCPs_and_constant_average_daily_CCPs(*ccp_data)
    plt.show()


if __name__ == "__main__":
    main()
