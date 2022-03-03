#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
from datetime import datetime, timezone
import csv
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


def get_SARSCoV2_peaks(seeds, CADS, WCADS):
    print('    ', 'seed', 'score', 'emax', 'cmax', 'smax',
          'emax_i', 'cmax_i', 'smax_i', 'emax_date', 'cmax_date', 'smax_date',)

    def get_data(dictionary, dict_name):
        scores = []
        emaxs = []
        cmaxs = []
        smaxs = []
        emaxs_index = []
        cmaxs_index = []
        smaxs_index = []
        emaxs_date = []
        cmaxs_date = []
        smaxs_date = []
        secmaxs = []
        for n, (seed, v) in enumerate(dictionary.items()):
            pdates = v['pdates']
            scores.append(v['score'][1])
            emp = np.array(v['plncases'], dtype=np.int16)
            covid19 = np.array(v['covid_days_ncases'], dtype=np.int16)
            sarscov2 = np.array(v['inf_days_ncases'], dtype=np.int16)
            emaxs.append(emp.max())
            cmaxs.append(covid19.max())
            smaxs.append(sarscov2.max())
            emaxs_index.append(np.where(emp == emaxs[n])[0][0])
            cmaxs_index.append(np.where(covid19 == cmaxs[n])[0][0])
            smaxs_index.append(np.where(sarscov2 == smaxs[n])[0][0])
            emaxs_date.append(pdates[emaxs_index[n]])
            cmaxs_date.append(pdates[cmaxs_index[n]])
            smaxs_date.append(pdates[smaxs_index[n]])
            # 1st May 2020 Secondary Peak
            sec_date = convert_str_to_datetime('01/05/2020')
            sec_index = pdates.index(sec_date)
            secmaxs.append(sarscov2[sec_index])
        for n, seed in enumerate(seeds):
            print(f'{dict_name}:', seed, scores[n],
                  emaxs[n], cmaxs[n], smaxs[n],
                  emaxs_index[n], cmaxs_index[n], smaxs_index[n],
                  emaxs_date[n].date(), cmaxs_date[n].date(),
                  smaxs_date[n].date())
        print('Secondary SARSCoV2 Peak: ', secmaxs)

    get_data(CADS, 'CADS')
    get_data(WCADS, 'WCADS')

    # print('CAD:', seed, score, emp.max(), covid19.max(), sarscov2.max())
    # for n, (seed, v) in enumerate(WCADS.items()):
    #     score = round(float(v['score'][1]), 2)
    #     emp = np.array(v['plncases'], dtype=np.int16)
    #     covid19 = np.array(v['covid_days_ncases'], dtype=np.int16)
    #     sarscov2 = np.array(v['inf_days_ncases'], dtype=np.int16)
    #     print('WCAD:', seed, score, emp.max(), covid19.max(), sarscov2.max())


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
    get_SARSCoV2_peaks(seeds, CADS, WCADS)
