import numpy as np
from numpy.random import Generator, PCG64DXSM, SeedSequence
from scipy.stats import norm

from pathlib import Path
from datetime import datetime, timezone
from itertools import count, repeat
from math import ceil
from operator import itemgetter
import statistics as stats
import csv
import concurrent.futures as cf
import threading
from debug_functions import prv, pri, npri


def get_estimated_mus_mean_median(*sources):
    # Read in predicted mus results to get x-axis and y-axis data
    predicted_mus = []
    xdates = []
    for source in sources:
        pmus, xd = read_predicted_mus(source)
        predicted_mus.append(pmus.values())
        xdates.append(xd)
    # xdates is x-axis data
    estimated_mus = []
    for per_day in zip(*predicted_mus):
        daily_samples = []
        for samples in per_day:
            daily_samples.extend(samples)
        estimated_mus.append(daily_samples)
    # Get the mean and median of the predicted mus
    ymus_mean = get_mean(estimated_mus)
    ymus_median = get_median(estimated_mus)
    return xdates[0], ymus_mean, ymus_median


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


def now():
    '''Function returns currnet date and time in "dd/mm/yyyy hh:mm:ss" format.
    '''
    return datetime.now().strftime("%d/%m/%Y %H:%M:%S")


def get_mean_and_median(pmus):
    ymus = list(pmus)
    mean = get_mean(ymus)
    median = get_median(ymus)
    return ymus, mean, median


def get_mean(y):
    '''Get the mean of y.'''
    ymean = []
    for i in y:
        try:
            mm = ceil(stats.mean(i))
        except stats.StatisticsError:
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
        except stats.StatisticsError:
            ymedian.append(0)
        else:
            ymedian.append(mm)
    return ymedian


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


class Predict_SARSCOV2_COVID19_VIA_CPTPP():
    '''Class to estimate the SARS-CoV-2 and COVID-19 trend using the Cartesian
    Product Transpose PP (CPTPP) method and CAD and WCAD to measure
    resemblance. Estimates with lowest CAD and WCAD values are updated to the
    CAD and WCAD .csv files.

    Inputs from previous backcasting-forecasting algorithm:
      xdates - list; dates of earliest possible SARS-CoV-2 infection to last
                     COVID-19 confirmation
      ymus - list; daily mean COVID-19 Confirmation Period predicted by
                   backcasting-forecasting(i.e. mus)

    Other inputs:
      covid19_csvfile - Path; the .csv file containing empical data of imported
                              and local COVID-19 cases
      mu_min_max_step - tuple; the range of mus (min, max, interval)
      sigma_constants - tuple; curve fitting constants for sigma as a function
                               of mus
      debug - boolean; Debugging default to False
    '''

    def __init__(self, xdates, ymus, covid19_csvfile, mu_min_max_step,
                 sigma_constants, debug=False, prefix='s'):
        if debug:
            # pri('xdates', xdates)
            # pri('ymus', ymus)
            prv('covid19_csvfile', covid19_csvfile)
            pri('mu_min_max_step', mu_min_max_step)
            pri('sigma_constants', sigma_constants)
            prv('debug', debug)
        self.debug = debug
        self.cwd = Path(__file__).parent
        epidemic_data = self._get_daily_cases(covid19_csvfile)
        self.imported = epidemic_data[0]
        self.local = epidemic_data[1]
        self.mu_sigma = self._get_mu_sigma(mu_min_max_step, sigma_constants)
        self.pdates = xdates
        self.pdays = self._get_local_days(xdates, self.local)
        self.plncases = self._get_local_ncases(xdates, self.local)
        self.pmus = ymus
        self.nmissing = 5
        self.cptpp_mus_arrays_list = \
            self._get_cptpp_mus_arrays_list(mu_min_max_step, self.nmissing)
        self.index0_top5, self.index0_others = \
            self._get_missing_pmus_ranked_index(self.pmus, self.plncases)
        self.cadfiles, self.wcadfiles = self.out_file_names(prefix)
        if debug:
            prv('self.debug', self.debug)
            pri('self.mu_sigma', self.mu_sigma)
            pri('self.pdates', self.pdates)
            pri('self.pdays', self.pdays)
            pri('self.plncases', self.plncases)
            pri('self.pmus', self.pmus)
            prv('self.nmissing', self.nmissing)
            pri('self.cptpp_mus_arrays_list', self.cptpp_mus_arrays_list)
            pri('self.index0_top5', self.index0_top5)
            pri('self.index0_others', self.index0_others)
            print('\n', f'self.cadfiles -> {len(self.cadfiles)}')
            for k, v in self.cadfiles.items():
                for i in v:
                    print(k, i)
            print('\n', f'self.wcadfiles -> {len(self.wcadfiles)}')
            for k, v in self.wcadfiles.items():
                for i in v:
                    print(k, i)
            # self._plot_pmus()
        pri('self.index0_top5', self.index0_top5)
        pri('self.index0_others', self.index0_others)

    def out_file_names(self, prefix):
        '''Returns two dictionaries {steps:[filenames]} for CAD and WCAD.'''
        cads = {}
        wcads = {}
        start = mu_min_max_step[0]
        stop = mu_min_max_step[1]+1
        interval = mu_min_max_step[2]
        for i in self.index0_top5.keys():
            if i == 0:
                cads[i] = [self.cwd/Path(f'{prefix}_stp{i}_o{mu}_cad.csv')
                           for mu in range(start, stop, interval)]
                wcads[i] = [self.cwd/Path(f'{prefix}_stp{i}_o{mu}_wcad.csv')
                            for mu in range(start, stop, interval)]
            else:
                cads[i] = [self.cwd/Path(f'{prefix}_stp{i}_cad_cad.csv'),
                           self.cwd/Path(f'{prefix}_stp{i}_wcad_cad.csv')]
                wcads[i] = [self.cwd/Path(f'{prefix}_stp{i}_cad_wcad.csv'),
                            self.cwd/Path(f'{prefix}_stp{i}_wcad_wcad.csv')]
        return cads, wcads

    def _get_daily_cases(self, csvfile):
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

    def _get_mu_sigma(self, mu_min_max_step, sigma_constants):
        '''Function to define the average number of days taken to confirm
        local covid19 cases and their corresponding standard deviations.

        Arguments:
         mu_min_max_step - min. and max. average number of days taken to
                           locally confirm COVID19 after SARS-CoV2 infection.
         sigma_constants - curve fitting constants (a,b,c) for a quadratic
                           mu-sigma relationship.

        Return:
         mu_sigma - a dictionary of the mu-sigma relationship
                    key = a range of the average number of days taken to
                          locally confirm COVID19 after SARS-CoV2 infection
                    value = standard deviations of the key
        '''
        mu_min, mu_max, _ = mu_min_max_step
        a, b, c = sigma_constants
        mu = [i for i in range(mu_min, mu_max+1, 1)]
        mu_sigma = {x: a*x**2 + b*x + c for x in mu}
        return mu_sigma

    def _get_local_days(self, xdates, local):
        '''Return a list of the day numbers where day 0 denotes the day that the
        1st local COVID19 case was confirmed. The len of this list is same as
        xdates

         args:
          xdates - a list of the dates of the predicted mus.
          local  - a dict; stores {date: no. of local cases}
        '''
        local_covid0day_date = next(
            (kv[0] for n, kv in enumerate(local.items()) if kv[1]), None)
        # prv('local_covid0day_date',local_covid0day_date)
        covid0day_index_in_xdates = next(
            (n for n, x in enumerate(xdates) if x == local_covid0day_date),
            None)
        # prv('covid0day_index_in_xdates', covid0day_index_in_xdates)
        xdays = list(range(-covid0day_index_in_xdates,
                           len(xdates)-covid0day_index_in_xdates
                           )
                     )
        # pri( 'xdays', xdays)
        # for n, i in enumerate(xdates):
        #    print(xdays[n], i)
        return xdays

    def _get_local_ncases(self, xdates, local):
        '''Return a list of the day numbers where day 0 denotes the day that the
        1st local COVID19 case was confirmed. The len of this list is same as
        xdates

         args:
          xdates - a list of the dates of the predicted mus.
          local  - a dict; stores {date: no. of local cases}
        '''
        local_1st_date = list(local.keys())[0]
        # prv('local_1st_date',local_1st_date)
        local_1st_date_index_in_xdates = next(
            (n for n, x in enumerate(xdates) if x == local_1st_date), None)
        # prv('local_1st_date_index_in_xdates', local_1st_date_index_in_xdates)
        ncases = [0 for _ in range(local_1st_date_index_in_xdates)] + \
            list(local.values())
        # pri( 'ncases', ncases)
        # for n, i in enumerate(xdates):
        #    print(ncases[n], i)
        return ncases

    def _get_cptpp_mus_arrays_list(self, mu_min_max_step, nmissing):
        # Get a list of Numpy array of mus
        musmin = mu_min_max_step[0]
        musmax = mu_min_max_step[1]
        step = mu_min_max_step[2]
        mus_array = np.array(list(range(musmin, musmax+1, step)),
                             dtype=np.uint64)
        mus_arrays_list = [mus_array] * nmissing
        # npri('mus_array', mus_array)
        # pri('mus_arrays_list', mus_arrays_list)
        return mus_arrays_list

    def _get_missing_pmus_ranked_index(self, pmus, plncases):
        # Get index of estimated mus that are missing mu value)
        index0 = [i for i, j in zip(count(), pmus) if j == 0]
        # pri('index0', index0)
        # Get the empirical no. of daily COVID19 cases for those days missing
        # mu value.
        L = [plncases[i] for i in index0]
        # pri('L', L)
        # Get the index of elements of L when L is sorted in decending order.
        index_of_ranked_L = [t[0] for t in sorted(enumerate(L),
                                                  key=itemgetter(1),
                                                  reverse=True)]
        # pri('index_of_ranked_L', index_of_ranked_L)
        # Get number of analysis step. Either it is floor divsion of the no. of
        # missing mu elements or 1 more if floor division has no remainder.
        steps = len(index0) // self.nmissing
        remainder = len(index0) % self.nmissing
        if remainder != 0:
            steps += 1
        # print(steps, remainder)
        # Set up array
        index0_top5 = {}
        index0_others = {}
        for n, ss in enumerate(range(0, len(index0), self.nmissing)):
            index_of_top5 = index_of_ranked_L[ss:ss+self.nmissing]
            index_of_others = index_of_ranked_L[ss+self.nmissing:]
            index0_top5[n] = [index0[i] for i in index_of_top5]
            index0_others[n] = [index0[i] for i in index_of_others]
            pri('index0_top5', index0_top5)
            pri('index0_others', index0_others)
        return index0_top5, index0_others

    def start1plus(self, seed=None, cpus=1, c=1, d=1):
        '''Method to analyse COVID19 confirm case data asynchronously to find
        the lowest CAD and WCAD of all possibile sequencing of mu that is
        generated by Cartesian Product for the missing elements of mus.

        kwargs:
          seed - an integer; NumPy random number generator seed.
          cpus - an integer; number of cpu to use.
          c, d - constants for Cartesian Product algorithm.

        Continuous Concurrent Job Submission Algorithm:
        To submit a continuous stream of jobs and analyse them concurrently,
        the Python's concurrent.futures.ProcessPoolExecutor is used within two
        while-loops. The function of the inner while-loop is to submit jobs
        continously to a pool of cpu cores. Their submission is regulated by
        the "job_submit" event flag; the total number of jobs submit in each
        loop can never exceed the maximum number of cpu core in the cpu core
         pool. The function of the outer while-loop is to ensure the resources
        shared between the "MainThread" (which is used to submit job) and a
        "Thread-1" (which is used to post-process the future of the submitted
        jobs) are not deadlocked. Also, this loop monitors the event flag
        "end_concurrent_execution". When its value is True, the outer
        while-loop will be exited and hence terminate the continuous streaming
        of concurrent job submissions.
        '''
        # 1. Localise class attributes
        debug = self.debug
        pdays = self.pdays
        pdates = self.pdates
        plncases = self.plncases
        mu_sigma = self.mu_sigma
        mu_range = list(self.mu_sigma.keys())
        steps = list(self.cadfiles.keys())

        # 2 Perform Resemblance Analysis for Step 1 onwards
        for step in steps[1:]:
            # 2.1 Check that all CAD and WCAD files in previous Step exists
            #     before performaning resemblance analysis.
            ps_cadexists = [file.exists() for file in self.cadfiles[step-1]]
            ps_wcadexists = [file.exists() for file in self.wcadfiles[step-1]]
            if False in ps_cadexists and False in ps_wcadexists:
                print(f'ps_cadexists={ps_cadexists}')
                print(f'ps_wcadexists={ps_wcadexists}')
                print(f'Step {step-1} is not performed or is incomplete.'
                      f'Please complete Step {step-1} first!')
                exit()
            else:
                print('All CAD and WCAD files exists!\n')
            # 2.2 Get index0 of top5 for current step
            top5 = self.index0_top5[step]
            print(f'step={step}')
            print(f'top5={top5}')
            # 2.3 Get index0 of top5 for current step
            cptpp_mus_arrays_list = [np.array(mu_range)] * len(top5)
            # 2.4 Get pmus of lowest CAD and lowest WCAD of previous step

            def get_pmus_of_least_cad_wcad_in_previous_step(
                    step_cadfiles, step_wcadfiles):
                min_pmus = []  # stores the pmus of the least CAD and WCAD
                step_files = (step_cadfiles, step_wcadfiles)
                for files in step_files:
                    if step == 1:
                        # Get all CAD/WCAD results of step
                        results = {mu_range[n]: read_CAD_WCAD(file)
                                   for n, file in enumerate(files)}
                        for k, v in results.items():
                            for ll, w in v.items():
                                print(k, '  ', ll)
                        # Get CAD and WCAD scores and final mus
                        scores = np.array([results[mu]['score'][1]
                                           for mu in mu_range])
                        final_mus = np.array([results[mu]['final_mus']
                                              for mu in mu_range])
                        # get Index of least CAD/WCAD scores
                        min_score = np.amin(scores)
                        min_score_index = (np.where(scores ==
                                                    min_score)[0].item())
                        # get predicted mus of least CAD/WCAD scores
                        min_pmus.append(final_mus[min_score_index])
                        if debug:
                            npri('scores', scores)
                            print('min_score', min_score)
                            print('min_score_index', min_score_index)
                    else:
                        # Get all CAD/WCAD results of step
                        results = {n: read_CAD_WCAD(file)
                                   for n, file in enumerate(files)}
                        for k, v in results.items():
                            for ll, w in v.items():
                                print(k, '  ', ll)
                        # Get CAD and WCAD final mus
                        final_mus = np.array([results[i]['final_mus']
                                              for i in range(2)])
                        # get predicted mus of least CAD/WCAD scores
                        min_pmus.append(final_mus)
                return min_pmus

            min_cad_wcad_pmus = get_pmus_of_least_cad_wcad_in_previous_step(
                self.cadfiles[step-1], self.wcadfiles[step-1])
            if debug:
                for pmus in min_cad_wcad_pmus:
                    pri('pmus', pmus)
            # 2.5 Find CAD and WCAD for least CAD and least WCADS
            # nR=0 is for least CAD of previous step
            # nR=1 is for least WCAD of previous step
            for nR, pmusR in enumerate(min_cad_wcad_pmus):
                if step == 1:
                    pmus = pmusR
                else:
                    pmus = pmusR[nR]
                # File to write CAD and WCAD results
                cadfile = self.cadfiles[step][nR]
                wcadfile = self.wcadfiles[step][nR]
                # concurrent.future result variables
                lcad = {'batch': 0, 'score': (0, 0), 'final_mus': None,
                        'inf_days_bins': None, 'inf_days_ncases': None,
                        'covid_days_ncases': None, 'batch_size': None,
                        'plncases': None, 'pdates': None, 'seed': None}
                lwcad = {'batch': 0, 'score': (0, 0), 'final_mus': None,
                         'inf_days_bins': None, 'inf_days_ncases': None,
                         'covid_days_ncases': None, 'batch_size': None,
                         'plncases': None, 'pdates': None, 'seed': None}

                def _init_seed():
                    '''Initize seed for NumPy Psuedo-Random Generator. If None,
                    current date_time is used.
                    '''
                    if not seed:  # use current date and time
                        return int(
                            str(datetime.now().timestamp()).replace('.', ''))
                    else:
                        return seed

                def _get_cartesian_product_transpose_pp_gen(a, c, d):
                    '''Returns a generator of Paul Panzer's Cartesian Product
                    Transpose PP function. The original
                    cartesian_product_transpose_pp function can't handle
                    problem size that is too large. Discretizing the function's
                    operation through two for-loops and using yielding to
                    create a generator can yield smaller chunks of the
                    cartesian product hence be able to handle larger problems.
                    args:
                      a - a list of NumPy arrays
                      c - the number of leading dimensions to chop up first
                      d - the number of leading dimensions to chop up second
                    '''
                    def cartesian_product_transpose_pp(arrays):
                        '''Returns the cartesian product of nested arrays.
                        args:
                          arrays - a list of NumPy arrays
                        '''
                        rows = len(arrays)
                        dtype = np.result_type(*arrays)
                        arr = np.empty((rows, *map(len, arrays)), dtype=dtype)
                        idx = slice(None), *repeat(None, rows)
                        for i, a in enumerate(arrays):
                            arr[i, ...] = a[idx[:rows-i]]
                        Tarr = arr.reshape(rows, -1).T
                        return Tarr

                    for x in cartesian_product_transpose_pp(a[:c]):
                        for y in cartesian_product_transpose_pp(
                                [*x[:, None], *a[c:d]]):
                            yield cartesian_product_transpose_pp(
                                [*y[:, None], *a[d:]])

                def _fix_missing_pmus_with_top5_only(top5_index, pmus,
                                                     batch_cptpp):
                    array2d = np.array([pmus]*len(batch_cptpp), dtype=np.int16)
                    batch_cptpp = np.array(batch_cptpp, dtype=np.int16)
                    array2d[:, top5_index] = batch_cptpp[:, ::-1]
                    if debug:
                        pri('top5_index', top5_index)
                        pri('pmus', pmus)
                        pri('batch_cptpp', batch_cptpp)
                        pri('after array2d', array2d)
                    return array2d

                def _update_lcad_lwcad(results, type=None):
                    nonlocal lcad, lwcad
                    results.append(pdates)
                    results.append(seed)
                    # pri('results', results)
                    if type in ['cad', 'wcad']:
                        if type in 'cad':
                            for n, k in enumerate(lcad.keys()):
                                lcad[k] = results[n]
                        elif type in 'wcad':
                            for n, k in enumerate(lwcad.keys()):
                                lwcad[k] = results[n]
                    else:
                        raise ValueError(f'type must be "cad" or "wcad". '
                                         f'{type} is invalid')

                def _write_results(results, target='out.csv'):
                    '''Function to write the progress of the prediction into
                     file.
                     '''
                    header = list(results.keys())
                    # print('header', header)
                    with open(target, mode='w', encoding='utf-8') as f:
                        csv_writer = csv.writer(f, quoting=csv.QUOTE_NONE)
                        csv_writer.writerow(header)  # write header
                        for h in header:
                            v = results[h]
                            # print(h, v)
                            if isinstance(v, int):
                                v = [v]
                            csv_writer.writerow(v)  # write each row

                def _get_job_fn_args(top5):
                    '''Function to get the arguments needed by function
                    self.analyse().
                    '''
                    nonlocal batch, batch_size, seedsequence
                    child_seed = seedsequence.spawn(1)
                    batch = next(batchcounter)
                    try:
                        batchcptppmus = next(cptppgen)
                    except StopIteration:
                        end_concurrent_execution.set()
                        print(
                            f'\nStopIteration : All '
                            f'{len(cptpp_mus_arrays_list)*len(mu_sigma)}'
                            f' possibilties considered.')
                        return False
                    else:
                        # print('batchcptppmus',batchcptppmus)
                        batchfixedpmus = _fix_missing_pmus_with_top5_only(
                            top5, pmus, batchcptppmus)
                    batch_size = len(batchfixedpmus)
                    batch_cs = [child_seed[0].entropy] * batch_size
                    # batch_rng = [default_rng(s) for s in batch_cs]
                    batch_rng = [Generator(PCG64DXSM(s)) for s in batch_cs]
                    if debug:
                        prv('batch', batch)
                        pri('batchfixedpmus', batchfixedpmus)
                        prv('batch_size', batch_size)
                        pri('batch_cs', batch_cs)
                        pri('batch_rng', batch_rng)
                        pri('pdays', pdays)
                        pri('plncases', plncases)
                        pri('mu_sigma', mu_sigma)
                        prv('cadfile', cadfile)
                        prv('wcadfile', wcadfile)
                    return (batch, batch_rng, batchfixedpmus, pdays, plncases,
                            mu_sigma, debug)

                def _process_future(result):
                    '''Function to process done future.

                    Local variable:
                     predicted_mus - a dict(); key -> day number
                                               value -> a list of predicted mus
                     xdates - a list of datetime.datetime() where ech item
                              stores the dates of the day number.
                    '''
                    nonlocal lcad, lwcad
                    # 1. Update and write lcad & lwcad results
                    cad_res = result[0]
                    wcad_res = result[1]
                    if cad_res[0] == 1 or wcad_res[0] == 1:
                        # print('\n1st result')
                        # 1st result
                        _update_lcad_lwcad(cad_res, type='cad')
                        _update_lcad_lwcad(wcad_res, type='wcad')
                        _write_results(lcad, target=cadfile)
                        _write_results(lwcad, target=wcadfile)
                        # print('1st result ', cad_res[0:2])
                        # print('1st result', wcad_res[0:2])
                    else:
                        # subsequent results:
                        # compare results and store results only if new
                        # least_score is lower print('\nsubsequent results')
                        if cad_res[1][1] < lcad['score'][1]:
                            _update_lcad_lwcad(cad_res, type='cad')
                            _write_results(lcad, target=cadfile)
                            print('## Updated CAD ', cad_res[0:2])
                        if wcad_res[1][1] < lwcad['score'][1]:
                            _update_lcad_lwcad(wcad_res, type='wcad')
                            _write_results(lwcad, target=wcadfile)
                            print('## Updated WCAD ', wcad_res[0:2])
                    max_possibilities = len(mu_sigma)**len(
                        cptpp_mus_arrays_list)
                    if cad_res[0]*cad_res[6] >= max_possibilities:
                        return True
                    if wcad_res[0]*wcad_res[6] >= max_possibilities:
                        return True

                def _on_complete(future):
                    '''Callback to process future when done.'''
                    nonlocal jobs_completed
                    # Process future
                    if _process_future(future.result()):
                        end_concurrent_execution.set()
                    # Discard completed future
                    running_futures.discard(future)
                    # Set 'job_submit' event flag to True when the number of
                    # futures is less than the number of cpus
                    if len(running_futures) < cpus:
                        submit_job.set()
                    # Update record on the number of completed jobs
                    jobs_completed += 1
                    # Wake up all threads waiting on this condition.
                    with job_cond:
                        job_cond.notify_all()

                # -------------------------------------------------------
                # Setup NumPy SeedSequence to produce repeatable pseudo-random
                # numbers across multiple processes.
                seedsequence = SeedSequence(_init_seed())
                # Setup batch counter
                batchcounter = count(start=0, step=1)
                batch = next(batchcounter)
                batch_size = None
                # Get_cartesian product transpose pp generator that will create
                #  values for the missing estimates of the estimated mus.
                cptppgen = _get_cartesian_product_transpose_pp_gen(
                    cptpp_mus_arrays_list, c, d)
                # Setup concurrent.futures variables
                running_futures = set()
                jobs_completed = 0
                job_cond = threading.Condition()
                # This threading event flag is initially false.
                end_concurrent_execution = threading.Event()
                submit_job = threading.Event()
                submit_job.set()
                # Notes:
                # - The main cpu core "MainThread" will submit jobs to
                #   different cpus cores.
                # - The inner while-loop ensures continuous job submission.
                # - A "submit_job" event flag is used to ensure the number of
                #   jobs submitted do not exceed the max. number of cpus core
                #   that is to be used.
                # - When a job is completed, the main cpu core "Thread-1" will
                #   post-process the job's future.
                # - job_cond.wait() and job_cond.notify_all() are used to exit
                #   the job submission loop and prevent a deadlock scenario
                #   between the activities of the "MainThread" and "Thread-1".
                # - When the event flag end_concurrent_execution value is True,
                #   the concurrent submission of jobs to the cpu core pool will
                #   cease. executor will await all running_futures to be done
                #   before shutting down itself.
                # Start concurrent execution
                print()
                print(f'** Concurrent analyses of least CAD & WCAD for '
                      f'Step={step} nR={nR} **')
                with cf.ProcessPoolExecutor(max_workers=cpus) as executor:
                    while True:  # job condition loop
                        while True:  # job submission loop
                            # Submit job
                            job_fn = self.analyse_cptpp
                            job_fn_args = _get_job_fn_args(top5)
                            if job_fn_args:
                                future = executor.submit(job_fn, *job_fn_args)
                            else:
                                break
                            future.add_done_callback(_on_complete)
                            running_futures.add(future)
                            # Show submitted jobs in terminal
                            stats = f'{now()}  cpus={cpus}  '\
                                f'rf={len(running_futures)}  step={step}  '\
                                f'nR={nR}  batch={batch}  size={batch_size}  '\
                                f'possibilities={batch_size*batch}  '
                            print(stats)
                            # Define event condition to pause job submission
                            # loop
                            if len(running_futures) >= cpus:
                                submit_job.clear()
                            # Allow job loop to proceed only when the event
                            # "submit_job" has a True value.
                            submit_job.wait()
                        # Outside the job submission loop, here we check the
                        # Wait to be notified
                        try:
                            with job_cond:
                                job_cond.wait()
                            # releases the underlying lock and then blocks
                            # until it is awakened by notify_all().
                            # call for the same condition variable in another
                            # thread.
                        except KeyboardInterrupt:
                            print(f'\n## DETECTED {len(running_futures)}'
                                  f' INCOMPLETED FUTURES ##')
                            print(' ...ATTEMPTING TO CANCEL THESE FUTURES')
                            print(' ...AFTER THESE FUTURES COMPLETES, PRESS '
                                  'CTRL+C MULTIPLE TIMES UNTIL SCRIPT ENDS.\n')
                            for future in running_futures:
                                _ = future.cancel()
                            # Ensure concurrent.futures.executor jobs finishes.
                            _ = cf.wait(running_futures, timeout=None)
                        # Condition to terminate the continuous streaming of
                        # job submissions.
                        if end_concurrent_execution.is_set():
                            print()
                            print(f'** Completed Concurrent analyses of least'
                                  f' CAD & WCAD for Step={step} nR={nR} **')
                            break
            print()
            print(f'** Get pmus of the lowest CAD & WCAD of Step {step}')

    def analyse_cptpp(self, batch, batchrng, batchmus, pdays, pncases,
                      mu_sigma, debug):
        '''Function to estimate the SARS-CoV-2 and Local COVID-19 trends using
        the daily mean covid19 confirmation periods, i.e. ??.
        '''
        if debug:
            print('\n\n #### analyse_cptpp ####\n')

        def _get_norm_rsv(mus, sigmas, size, rg):
            castdays = np.int16(
                norm.rvs(loc=mus, scale=sigmas, size=size, random_state=rg)
            )
            return castdays

        def _fix_predictions_ends(ptype, projections, current,
                                  samples_mu, samples_sigma, rng,
                                  imported0day, locallastday):
            projections_mins = np.min(projections, axis=1)
            projections_maxs = np.max(projections, axis=1)
            if debug:
                npri('projections_mins', projections_mins)
                npri('projections_maxs', projections_maxs)

            exceeds = np.where(
                (projections < imported0day) | (projections > locallastday)
            )
            if debug:
                pri('exceeds', exceeds)
                prv('exceeds[0].size', exceeds[0].size)
            if exceeds[0].size > 0:
                if debug:
                    print('if exceeds is true')
                for r, c in zip(*exceeds):
                    while True:
                        # print(r, c, samples_mu[r,c], samples_sigma[r,c],
                        #      current[r,c])
                        cast_day = _get_norm_rsv(samples_mu[r, c],
                                                 samples_sigma[r, c],
                                                 1, rng[r])
                        if ptype in 'backcast':
                            day = current[r, c] - cast_day
                        elif ptype in 'forecast':
                            day = current[r, c] + cast_day
                        # print('cast_day, day, imported0day locallastday',
                        #      cast_day, day, imported0day, locallastday)
                        if day >= imported0day and day <= locallastday:
                            projections[r, c] = day[0]
                            break
            projections_mins = np.min(projections, axis=1)
            projections_maxs = np.max(projections, axis=1)
            if debug:
                npri('projections_mins', projections_mins)
                npri('projections_maxs', projections_maxs)
            return projections

        def _get_histogram(a, bins=10, range=None, normed=None, weights=None,
                           density=None):
            '''Returns the values of a histogram.'''
            hist, _ = np.histogram(a, bins=bins, range=range, normed=normed,
                                   weights=weights, density=density)
            return hist

        def _get_lowest_cad_and_wcad(lncases, encases):
            '''Get location of the lowest Cumulative Absolute Difference (CAD)
            and Weighted Cumulative Absolute Difference (WCAD) of the estimated
            daily local COVID19 case trends.

            args:
              lncases - numpy.array(); empirical daily number of local covid19
                                       cases
              encases - numpy.array(); estimated daily number of local covid19
                                       cases

            return:
              A 2 element tuple. Each element is also a 2 element tuple.
              1st element - (batch row with lowest AD, lowest CAD)
              2nd element - (batch row with lowest WAD, lowest WCAD)
            '''
            ad = np.absolute(encases - lncases)
            ad_sum = np.sum(ad, axis=1)
            ad_sum_min = ad_sum.min()
            if debug:
                npri('ad', ad)
                npri('ad_sum', ad_sum)
                prv('ad_sum_min', ad_sum_min)

            weights = lncases / lncases.sum()
            wad = weights * ad
            wad_sum = np.sum(wad, axis=1)
            wad_sum_min = wad_sum.min()
            if debug:
                npri('weights', weights)
                npri('wad', wad)
                npri('wad_sum', wad_sum)
                prv('wad_sum_min', wad_sum_min)

            ad_sum_min_index = np.where(ad_sum == ad_sum_min)
            wad_sum_min_index = np.where(wad_sum == wad_sum_min)
            if debug:
                pri('ad_sum_min_index', ad_sum_min_index)
                pri('wad_sum_min_index', wad_sum_min_index)

            return ((int(ad_sum_min_index[0][0]), ad_sum_min),
                    (int(wad_sum_min_index[0][0]), wad_sum_min))

        # -- start of analyse_cptpp ---
        covid_earliest_day = pdays[0]
        locallastday = pdays[-1]
        batch_size = len(batchmus)
        axis0 = batch_size
        axis1 = len(pncases)
        saxis1 = np.sum(pncases).tolist()
        pshape = (axis0, axis1)
        if debug:
            prv('batch', batch)
            prv('covid_earliest_day', covid_earliest_day)
            prv('locallastday', locallastday)
            prv('batch_size', batch_size)
            pri('batchrng', batchrng)
            prv('axis0', axis0)
            prv('axis1', axis1)
            pri('pshape', pshape)

        # Create final 2D arrays
        final_mus = np.array(batchmus, dtype=np.int64)
        final_sigmas = np.array([[mu_sigma[j] for j in i] for i in final_mus])
        final_days = np.broadcast_to(pdays, pshape)
        if debug:
            npri('final_mus', final_mus)
            npri('final_sigmas', final_sigmas)
            pri('pdays', pdays)
            npri('final_days', final_days)

        # Convert daily mu and sigmas and days for every COVID19 case
        samples_mu = np.repeat(final_mus, pncases, axis=1)
        samples_sigma = np.repeat(final_sigmas, pncases, axis=1)
        samples_days = np.repeat(final_days, pncases, axis=1)
        if debug:
            prv('sum(pncases)', sum(pncases))
            npri('samples_mu', samples_mu)
            npri('samples_sigma', samples_sigma)
            npri('samples_days', samples_days)

        # Backcast: predict when SARS-Cov2_Infections occurred for various mus
        #           and sigmas
        bcd = []
        for n, rng in enumerate(batchrng):
            bcd.append(
                _get_norm_rsv(samples_mu[n, :], samples_sigma[n, :], saxis1,
                              rng)
            )
        backcast_days = np.vstack(bcd).reshape(axis0, saxis1)
        inf_days = samples_days - backcast_days
        if debug:
            pri('bcd', bcd)
            npri('backcast_days', backcast_days)
            npri('inf_days', inf_days)
        inf_days = _fix_predictions_ends('backcast', inf_days, samples_days,
                                         samples_mu, samples_sigma, batchrng,
                                         covid_earliest_day, locallastday)
        if debug:
            npri('inf_days', inf_days)

        # Get daily number of predicted SARS-Cov2 infections
        inf_days_bins = np.arange(covid_earliest_day, locallastday+2, 1,
                                  dtype=np.int16)
        inf_days_ncases = np.apply_along_axis(_get_histogram, 1, inf_days,
                                              bins=inf_days_bins)
        itotaldays = np.sum(inf_days_ncases, axis=1)
        if debug:
            npri('inf_days_bins', inf_days_bins)
            npri('inf_days_ncases', inf_days_ncases)
            npri('itotaldays', itotaldays)

        # Forecast: predict when COVID19 was confirmed for various mus and
        #           sigmas
        fcd = []
        for n, rng in enumerate(batchrng):
            fcd.append(
                _get_norm_rsv(samples_mu[n, :], samples_sigma[n, :], saxis1,
                              rng)
            )
        forecast_days = np.vstack(fcd).reshape(axis0, saxis1)
        covid_days = inf_days + forecast_days
        if debug:
            pri('fcd', fcd)
            npri('forecast_days', forecast_days)
            npri('covid_days', covid_days)
        covid_days = _fix_predictions_ends('forecast', covid_days, inf_days,
                                           samples_mu, samples_sigma, batchrng,
                                           covid_earliest_day, locallastday)
        if debug:
            npri('covid_days', covid_days)

        # Get daily number of predicted COVID19 confirmed cases
        covid_days_bins = np.arange(covid_earliest_day, locallastday+2, 1,
                                    dtype=np.int16)
        covid_days_ncases = np.apply_along_axis(_get_histogram, 1, covid_days,
                                                bins=covid_days_bins)
        ctotaldays = np.sum(covid_days_ncases, axis=1)
        if debug:
            npri('covid_days_bins', covid_days_bins)
            npri('covid_days_ncases', covid_days_ncases)
            npri('ctotaldays', ctotaldays)

        # identify which row gives estimated COVID19 trend with lowest CAD and
        # WCAD
        empirical = np.array(pncases, dtype=np.int32)
        least_cad, least_wcad = _get_lowest_cad_and_wcad(empirical,
                                                         covid_days_ncases)
        if debug:
            npri('empirical', empirical)
            pri('least_cad', least_cad)
            pri('least_wcad', least_wcad)
        adrow = least_cad[0]
        wadrow = least_wcad[0]
        least_ad_res = [batch,
                        least_cad,
                        final_mus[adrow].tolist(),
                        inf_days_bins.tolist(),
                        inf_days_ncases[adrow].tolist(),
                        covid_days_ncases[adrow].tolist(),
                        batch_size,
                        pncases,
                        ]
        least_wad_res = [batch,
                         least_wcad,
                         final_mus[wadrow].tolist(),
                         inf_days_bins.tolist(),
                         inf_days_ncases[wadrow].tolist(),
                         covid_days_ncases[wadrow].tolist(),
                         batch_size,
                         pncases,
                         ]
        if debug:
            print('\nleast_ad_res = ', least_ad_res)
            print('\nleast_wad_res = ', least_wad_res)
        return least_ad_res, least_wad_res


if __name__ == "__main__":
    repo_dir = Path(__file__).resolve().parents[2]
    # Statistical ?? estimates of random seeds 1, 2, 3
    estimates_dir = repo_dir / '2_Results' / '1_statistical_??'
    e1 = estimates_dir / 'nd_predicted_mus_a1_s91023483923_pc8458.csv'
    e2 = estimates_dir / 'nd_predicted_mus_a2_s66576402094_pc8458.csv'
    e3 = estimates_dir / 'nd_predicted_mus_a3_s343090589475_pc8458.csv'
    xdates, ymus_mean, _ = get_estimated_mus_mean_median(e1, e2, e3)
    print(ymus_mean)
    # Inputs needed to predict SARS-CoV-2 and COVID-19 trends
    covid19_csvfile = 'COVID19_epidemic_trends.csv'
    mu_min_max_step = (1, 18, 1)
    sigma_constants = (-0.008665, 0.483888, 0.0)
    seed1 = 91023483923
    seed2 = 66576402094
    seed3 = 343090589475
    seeds = [seed1, seed2, seed3]
    labels = ['s1', 's2', 's3']
    cpus = 2
    c = 2    # Constant for cartesian_product_transpose_pp
    d = 3    # Constant for cartesian_product_transpose_pp
    # Define covid19_csvfile path
    covid_data = repo_dir / '4_Empirical_Data' / covid19_csvfile
    # Perform Step1+ of Modified Resemblance Algorithm for 3 random seeds
    for seed, label in zip(seeds, labels):
        analysis = Predict_SARSCOV2_COVID19_VIA_CPTPP(
            xdates, ymus_mean, covid_data, mu_min_max_step, sigma_constants,
            prefix=label
            # debug=True,
        )
        analysis.start1plus(seed=seed, cpus=cpus, c=c, d=d)
