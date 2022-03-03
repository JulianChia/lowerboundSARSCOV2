
from pathlib import Path
from datetime import datetime, timezone

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def read_sta_file(stafile):
    with open(stafile, mode='r') as file:
        lines = file.readlines()
        dt = []
        batch = []
        total = []
        iterations = []
        empty = []
        empty_pc = []
        incomplete = []
        incomplete_pc = []
        completed = []
        completed_pc = []
        for line in lines[1:]:
            a = line.split()
            dt.append(convert_str_to_datetime(a[0]+' '+a[1],
                                              dateformat="%d/%m/%Y %H:%M:%S",
                                              tzinfo=timezone.utc))
            batch.append(int(a[2]))
            iterations.append(int(a[3]))
            total.append(int(a[4]))
            empty.append(int(a[5]))
            empty_pc.append(int(a[5]) / int(a[4]) * 100)
            incomplete.append(int(a[7]))
            incomplete_pc.append(int(a[7]) / int(a[4]) * 100)
            completed.append(int(a[9]))
            completed_pc.append(int(a[9]) / int(a[4]) * 100)
    return (dt, batch, iterations, total, empty, empty_pc, incomplete,
            incomplete_pc, completed, completed_pc)


def convert_str_to_datetime(date, dateformat="%d/%m/%Y", tzinfo=timezone.utc):
    refdate = datetime.strptime(date, dateformat)
    refdate = refdate.replace(tzinfo=tzinfo)
    return refdate


def plot_performance(fname1, fname2, fname3):
    fontsize_axis = 16
    fontsize_axislabel = 20

    # Read in data from .sta files
    nd1 = read_sta_file(fname1)
    nd2 = read_sta_file(fname2)
    nd3 = read_sta_file(fname3)
    dt1, batch1, iterations1, total1, empty1, empty1_pc, incomplete1, incomplete1_pc, completed1, completed1_pc = nd1
    dt2, batch2, iterations2, total2, empty2, empty2_pc, incomplete2, incomplete2_pc, completed2, completed2_pc = nd2
    dt3, batch3, iterations3, total3, empty3, empty3_pc, incomplete3, incomplete3_pc, completed3, completed3_pc = nd3
    print('max iterations : ',
          max(iterations1), max(iterations2), max(iterations3),
          )
    print('max completed : ',
          max(completed1_pc), max(completed2_pc), max(completed3_pc),
          )
    print('min incomplete : ',
          min(incomplete1_pc), min(incomplete2_pc), min(incomplete3_pc),
          )
    print('min empty : ',
          min(empty1_pc), min(empty2_pc), min(empty3_pc),
          )

    # Start Time
    dt1_min = dt1[0]
    dt2_min = dt2[0]
    dt3_min = dt3[0]

    # Time difference in minutes w.r.t Start Time
    dt1_delta = [(i-dt1_min).total_seconds()/60 for i in dt1]
    dt2_delta = [(i-dt2_min).total_seconds()/60 for i in dt2]
    dt3_delta = [(i-dt3_min).total_seconds()/60 for i in dt3]
    print('max time : ', max(dt1_delta), max(dt2_delta), max(dt3_delta))

    # Construct plot
    fig, ax = plt.subplots(1, 2, sharex=False, sharey=False)

    # kwargs of plot lines
    alpha = 0.9
    col2 = 'coral'
    col3 = 'seagreen'
    col1 = 'deepskyblue'
    ls1 = 'solid'
    ls2 = 'dashed'
    ls3 = 'dotted'
    c1 = dict(label='ND1 Completed',   color=col1, linestyle=ls1, alpha=alpha)
    ic1 = dict(label='ND1 Incomplete', color=col2, linestyle=ls2, alpha=alpha)
    e1 = dict(label='ND1 None',        color=col3, linestyle=ls3, alpha=alpha)
    c2 = dict(label='ND2 Completed',   color=col1, linestyle=ls1, alpha=alpha)
    ic2 = dict(label='ND2 Incomplete', color=col2, linestyle=ls2, alpha=alpha)
    e2 = dict(label='ND2 None',        color=col3, linestyle=ls3, alpha=alpha)
    c3 = dict(label='ND3 Completed',   color=col1, linestyle=ls1, alpha=alpha)
    ic3 = dict(label='ND3 Incomplete', color=col2, linestyle=ls2, alpha=alpha)
    e3 = dict(label='ND3 None',        color=col3, linestyle=ls3, alpha=alpha)

    x_min_tup = (dt1_delta, dt1_delta, dt1_delta,
                 dt2_delta, dt2_delta, dt2_delta,
                 dt3_delta, dt3_delta, dt3_delta,)
    x_iter_tup = (iterations1, iterations1, iterations1,
                  iterations2, iterations2, iterations2,
                  iterations3, iterations3, iterations3,)
    y_pc_tup = (completed1_pc, incomplete1_pc, empty1_pc,
                completed2_pc, incomplete2_pc, empty2_pc,
                completed3_pc, incomplete3_pc, empty3_pc,)
    p_tup = (c1, ic1, e1,
             c2, ic2, e2,
             c3, ic3, e3,)

    # Plot Percentage vs Minutes
    for x0, y0, p in zip(x_min_tup, y_pc_tup, p_tup):
        ax[0].semilogx(x0, y0, **p)

    # Plot Percentage vs Iterations
    for x1, y1, p in zip(x_iter_tup, y_pc_tup, p_tup):
        ax[1].semilogx(x1, y1, **p)

    # Configure X-axes Label
    ax[0].get_yaxis().set_label_text(label='Percentage (%)', fontweight='bold',
                                     fontsize=fontsize_axislabel)
    ax[1].get_yaxis().set_label_text(label='Percentage (%)', fontweight='bold',
                                     fontsize=fontsize_axislabel)

    # Configure Y-axes Label
    ax[0].get_xaxis().set_label_text(label='Minutes', fontweight='bold',
                                     fontsize=fontsize_axislabel)
    ax[1].get_xaxis().set_label_text(label='Iterations', fontweight='bold',
                                     fontsize=fontsize_axislabel,)

    # Configure Scale of X-axis
    ax[0].set_ylim(0, 100)
    ax[1].set_ylim(0, 100)

    # Draw Ticks
    ax[0].tick_params(axis='both', which='both', labelsize=fontsize_axis)
    ax[1].tick_params(axis='both', which='both', labelsize=fontsize_axis)
    ax[0].yaxis.set_major_locator(ticker.MultipleLocator(10))
    ax[0].yaxis.set_minor_locator(ticker.MultipleLocator(1))
    ax[1].yaxis.set_major_locator(ticker.MultipleLocator(10))
    ax[1].yaxis.set_minor_locator(ticker.MultipleLocator(1))

    # Draw Grid
    ax[0].grid(linestyle='-', color='grey', which='major')
    ax[0].grid(linestyle=':', color='grey', which='minor', axis='x')
    ax[1].grid(linestyle='-', color='grey', which='major')
    ax[1].grid(linestyle=':', color='grey', which='minor', axis='x')

    # Draw Legend
    handles0, _ = ax[0].get_legend_handles_labels()
    handles1, _ = ax[1].get_legend_handles_labels()
    handles0 = handles0[0:3]
    handles1 = handles1[0:3]
    labels = ['Completed   - S1, S2, S3',
              'Incomplete   - S1, S2, S3',
              'No Estimate - S1, S2, S3']
    ax[0].legend(handles0, labels, loc='center right', fontsize=fontsize_axislabel)
    ax[1].legend(handles1, labels, loc='center right', fontsize=fontsize_axislabel)

    # Configure Subplot layout
    top = 0.980  # the top of the subplots of the figure
    bottom = 0.085  # the bottom of the subplots of the figure
    left = 0.050  # the left side of the subplots of the figure
    right = 0.980  # the right side of the subplots of the figure
    hspace = 0.150  # the amount of height reserved for space between subplots,
    # expressed as a fraction of the average axis height
    wspace = 0.130  # the amount of width reserved for space between subplots,
    # expressed as a fraction of the average axis width
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top,
                        wspace=wspace, hspace=hspace)

    # Draw Remarks
    msg0 = 'Completed'
    msg1 = 'Incomplete'
    msg2 = 'No Estimate'
    msg0x = 1.*1e2
    msg1x = 2.*1e1
    msg2x = 2
    msg0y = 83
    msg1y = 30
    msg2y = 10
    alpha = 0.6
    ax[0].text(msg0x, msg0y, msg0, ha="left", va="bottom", color='black',
               fontsize=16, alpha=alpha)
    ax[0].text(msg1x, msg1y, msg1, ha="left", va="bottom", color='black',
               fontsize=16, alpha=alpha)
    ax[0].text(msg2x, msg2y, msg2, ha="left", va="bottom", color='black',
               fontsize=16, alpha=alpha)
    msg0x = 8*1e4
    msg1x = 3*1e4
    msg2x = 8*1e2
    msg0y = 83
    msg1y = 30
    msg2y = 10
    ax[1].text(msg0x, msg0y, msg0, ha="left", va="bottom", color='black',
               fontsize=16, alpha=alpha)
    ax[1].text(msg1x, msg1y, msg1, ha="left", va="bottom", color='black',
               fontsize=16, alpha=alpha)
    ax[1].text(msg2x, msg2y, msg2, ha="left", va="bottom", color='black',
               fontsize=16, alpha=alpha)

    return fig


if __name__ == "__main__":
    # Get Data
    repo_dir = Path(__file__).resolve().parents[2]
    results_data = repo_dir / '2_Results' / '1_statistical_μ'
    fname1 = results_data / 'nd_predicted_mus_a1_s91023483923_pc8458.sta'
    fname2 = results_data / 'nd_predicted_mus_a2_s66576402094_pc8458.sta'
    fname3 = results_data / 'nd_predicted_mus_a3_s343090589475_pc8458.sta'
    # Show Data
    plot_performance(fname1, fname2, fname3)
    plt.show()
