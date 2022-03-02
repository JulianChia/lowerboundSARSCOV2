#!/usr/bin/env python3

"""
What does this module do?

1. Read in SG COVID-19 confirmed case population data (a .csv file).
2. Read in assumed quadratic mu-sigma relationship.
3. Use these info and proposed statistical algorithum to predict the probable
   SG SARS_CoV2 infections population trend that preceeded the SG COVID-19
   confirmed case population trend.


How to use this module?

1. You need to supply these info in the main function:
   -------------------
   Mandatory Arguments
   -------------------
    covid19_csvfile - str; the path of the SG COVID-19 confirmed case
                           population data .csv file
    mu_min - integer; the lowest average number of days to possibly confirm
                      COVID19 after SARS-CoV2 infection.
    mu_max - integer; the highest average number of days to possibly confirm
                      COVID19 after SARS-CoV2 infection.
    sigma_a and sigma_b - integers; constants to define quadratic relationship
                                    between the average number of days to
                                    confirm COVID19 (mu) and its standard
                                    deviation (sigma) in a normal distribution.
    earliest_arrival_date_of_imported_cases - str; yyyy-mm-dd format
   -------------------
    Optional Arguments
   -------------------
    seed - interger; to seed the random number generator of NumPy. Defaults to
                     current date-time.
    iterations - integer; the number of iterations analysed in each batch
                          analysis.
    cpus - interger; the number of cpu logical cores to analyse the data
                        concurrently. Defaults to all cpu cores.

2. The format of the covid19_csvfile .csv data file must be:
   1st row (header)           : Date, No. of Daily Imported cases,
                                No. of Daily Local cases
   2nd & subsequent rows(data): 2020-01-23, 2, 10
   Note: The format of the date in the first column must be yyyy-mm-dd.

3. This module will run the function initialize() and return a
   "Predict_COVID19_Confirmation_Periods" class object. Thereafter, it will run
   the "start()" method of that class object and write the predicted probable
   number of days to confirm COVID19 in an output file (default filename is
   "out.csv".)
"""
import numpy as np
from numpy.random import Generator, PCG64DXSM, SeedSequence
from scipy.stats import norm

from pathlib import Path
from itertools import count
from datetime import datetime, timedelta, timezone
import concurrent.futures as cf
import threading
import csv
import os

from debug_functions import prv, pri, npri


class Local_COVID19_Parameters():
    '''       Local_COVID19_Parameters( local_covid_dict )

    Class to store the vital data on confirmed Local COVID19 cases".

    Argument:
      local_covid_dict - Dictionary on confirmed Local COVID19 cases
                         {dates:number of cases}
    Attributes:
      self.vitals - a tuple; (first_date, last_date, period)
      self.dates  - a list; a range of COVID19 dates [1st case to last case]
      self.days   - a np.array; an enumeration of self.dates
      self.ncases - a np.array; the number of daily COVID19 cases in self.days
    '''

    def __init__(self, local_covid_dict):
        data = self.get_local_covid_data(local_covid_dict)
        self.vitals = data[0]
        self.dates = data[1]
        self.days = np.array(data[2])
        self.ncases = np.array(data[3])
        # debug
        pri('self.vitals', self.vitals)
        pri('self.dates', self.dates)
        npri('self.days',  self.days)
        npri('self.ncases', self.ncases)
        prv('self.ncases.cumsum()', self.ncases.cumsum())

    def get_local_covid_data(self, local_covid_dict):
        ''' Function to return key data on local_covid_dict. '''
        # Get 1st and last date and period of confirmed Local COVID19 cases
        first_date, first_date_num = self.find_first_case(local_covid_dict)
        last_date = list(local_covid_dict.keys())[-1]
        period = (last_date - first_date).days + 1
        # Define Local COVID19 cases date
        dates = list(local_covid_dict.keys())[first_date_num:]
        # Serialise the Local COVID19 cases date.
        # 0 denotes day in which the 1st case was reported.
        serial = [i for i in range(period)]
        # Get an array of the number of daily local COVID19 cases
        ncases = list(local_covid_dict.values())[first_date_num:]
        return (first_date, last_date, period), dates, serial, ncases

    def find_first_case(self, cases_dict):
        count = 0
        for k, v in cases_dict.items():
            if v > 0:
                first_case_date = k
                break
            count += 1
        return first_case_date, count


class Predict_COVID19_Confirmation_Periods():
    '''Predict_COVID19_Confirmation_Periods( mu_sigma, local,
    earliest_covid_date_possible=None )

    A class to predict the average number of days taken to confirm daily counts
    of Local COVID19 cases after the contraction of SARS-CoV2.

    Arguments:
      mu_sigma - a dict(); data on the assumed relationship between the average
                           number of days taken to confirm daily counts of
                           Local COVID19 cases (mu) and its standard deviation
                           (sigma).
      local    - an instance of Local_COVID19_Parameters().
      earliest_covid_date_possible - a datetime.datetime();
                                     earliest date possible when COVID19 could
                                     have been confirmed in Singapore.
    Attributes:
      self.mu_scenarios    - 1D np.array(); mu from mu_sigma
      self.sigma_scenarios - 1D np.array(); sigma from mu_sigma inserted with
                                            0.0 in the front.
      self.local           - see Arguments
      self.seed            - see Arguments
      self.covid_earliest_day
                           - an integer; day number of first reported imported
                                         COVID19 case.
      self.mus             - 2D np.array; for each mu_scenarios and day, stores
                                          the mu_scenarios element.
      self.sigmas          - 2D np.array; for each mu_scenarios and day, stores
                                          the corresponding standard deviations
                                          of self.mus.
    '''

    def __init__(self, mu_sigma, local, earliest_covid_date_possible=None):

        def _initialise_daily_mu(mu_scenarios, days):
            nrow = mu_scenarios.size
            ncol = days.size
            return np.full((nrow, ncol), mu_scenarios[:, None], dtype=np.uint8)

        def _initialise_daily_sigma(mu):
            '''Function returns the sigma corresponding to a given mu scenario.
            '''
            sigma = np.insert(self.sigma_scenarios[::-1], 0, 0.0)
            return sigma[mu]

        def _get_prediction_dates():
            '''Function returns the dates from firstday to lastday.'''
            numdays = local.days[-1].item() - self.covid_earliest_day.item()
            prv('numdays', numdays)
            date_list = [earliest_covid_date_possible + timedelta(days=x)
                         for x in range(numdays+1)]
            return date_list

        # Define the attributes to store read-in experimental data
        self.local = local
        if earliest_covid_date_possible:
            self.covid_earliest_day = local.days[0] - \
                (local.dates[0] - earliest_covid_date_possible).days
        else:
            self.covid_earliest_day = local.days[0]
        prv(' self.covid_earliest_day',  self.covid_earliest_day)

        # Define the attributes to store Mu and Sigma scenarios
        self.mu_scenarios = np.array(list(mu_sigma.keys()),
                                     dtype=np.uint8)[::-1]
        self.sigma_scenarios = np.array(list(mu_sigma.values()))[::-1]
        npri('self.mu_scenarios', self.mu_scenarios)
        npri('self.sigma_scenarios', self.sigma_scenarios)

        # Define the attributes to store the Daily Mu and Sigma 2D arrays
        self.mus = _initialise_daily_mu(self.mu_scenarios, local.days)
        self.sigmas = _initialise_daily_sigma(self.mus)
        npri('self.mus', self.mus)
        npri('self.sigmas', self.sigmas)

        # Define the attribute to store dates of covid19 analyse
        self.prediction_dates = _get_prediction_dates()
        pri('self.prediction_dates', self.prediction_dates)

    def start(self, seed=None, cpus=1, iterations=100, quota=500,
              out='predicted_mus', debug=False, factor=1.0,
              completion_target=100.0):
        '''Method to start analysing COVID19 confirm case data asynchronously.

        kwargs:
          seed - integer; NumPy random number generator seed.
          cpus - integer; number of cpu to use.
          iterations - an integer; number of iterations analysed by each cpu
                                   core or job batch.
          quota - integer; number of predictions to store for each day.
          out - str; file to store predictions.
          debug - bolean; True or False
          factor - float; accuracy factor. 1.0 denote predicted mu must be
                          identical to COVID19 confirm case data. 0.8 denotes
                          predicted mu must be within the range of being 80%
                          to 100% of COVID19 confirm case data.
          completion_target - float; 0.0 to 100.0

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

        def _init_debug():
            '''Function to ensure debug value is either True or False.'''
            if debug is True or debug is False:
                return debug
            else:
                raise ValueError('"debug" value must be either "True" or "False"')

        def _init_seed():
            '''Function to initize seed for NumPy Psuedo-Random Generator.
            If None, current date_time is used.'''
            if not seed:  # use current date and time
                return int(str(datetime.now().timestamp()).replace('.', ''))
            else:
                return seed

        def _get_job_fn_args():
            '''Function to get the arguments needed by function self.analyse()'''
            nonlocal batch, seedsequence
            child_seed = seedsequence.spawn(1)
            #print( 'child_seed = ', child_seed)
            batch = next(batchcounter)
            return (iterations, child_seed, batch, self.local.days,
                    self.local.ncases, self.mu_scenarios, self.sigma_scenarios,
                    self.covid_earliest_day, factor, quota)

        def _process_future(result):
            '''Function to process done future.

            Local variable:
             predicted_mus - a dict(); key -> day number
                                       value -> a list of predicted mus
             xdates - a list of datetime.datetime() where ech item stores the
                      dates of the day number.
            '''
            batch = result[1]
            iterations = result[2]
            # 1. Define predicted_mus and corresponding dates
            # if statusfile.stat().st_size <= 95: # for processing first future
            if statusfile.stat().st_size <= 191:  # for processing first future
                predicted_mus = result[0]
                xdates = self.prediction_dates
            else:
                predicted_mus, xdates = read_predicted_mus(source=resultfile)
                # Append future.result into existing predicted_mus.
                for k, v in result[0].items():
                    if len(predicted_mus[k]) < quota:
                        predicted_mus[k].extend(v)
            # 2. Write predictions to .csv file
            days = list(predicted_mus.keys())
            dates = [i.date() for i in xdates]
            mus_list = list(predicted_mus.values())
            headers = ['Days', 'Dates',
                       'Predicted Average Number of Days to Confirm COVID19 '
                       'After SARS-CoV2 Contraction']
            write_predicted_mus(days, dates, mus_list, headers,
                                target=resultfile)
            # FOR DEBUGGING: See status of predicted mus
            #mus_list_elemsize = [ len(i) for i in mus_list ]
            #index0 = [ i for i, j in zip(count(), mus_list) if len(j) == 0 ]
            # print( 'Number of Predictions Collected:\n',
            #       mus_list_elemsize, len(mus_list),
            #       '\nElements without a prediction:\n',
            #       mus_list_elemsize.count( 0 ), index0 )
            # 3. Return status of analysis
            prediction_status = [True if len(i) >= quota else
                                 False for i in mus_list]
            # True denotes completed
            # False denotes incomplete
            total = len(mus_list)
            empty = mus_list.count([])
            incomplete = prediction_status.count(False)
            completed = prediction_status.count(True)
            pc_empty = empty/total*100
            pc_incomplete = incomplete/total*100
            pc_completed = completed/total*100
            iterations = batch * iterations
            status = f'{now()}  {batch:6}  {iterations:>10}  {total:>6}  '\
                     f'{empty:>6}  {pc_empty:>.2f}%  '\
                     f'{incomplete:>6}  {pc_incomplete:>.2f}%  '\
                     f'{completed:>6}  {pc_completed:>.2f}%'
            write_status(status, target=statusfile)
            if False in prediction_status:
                # Prediction is completed
                return (True, pc_completed)
            else:
                # Prediction is not completed
                return (False, pc_completed)

        def _on_complete(future):
            '''Callback to process future when done.'''
            nonlocal jobs_completed, percentage_completed
            # Process future and set 'end_concurrent_execution' event flag to True
            # when sufficient data is obtained.
            quota_met, percentage_completed = _process_future(future.result())
            if quota_met:
                end_concurrent_execution.set()
            # Discard completed future
            running_futures.discard(future)
            # Set 'job_submit' event flag to True when the number of futures is
            # less than the number of cpus
            if len(running_futures) < cpus:
                submit_job.set()
            # Update record on the number of completed jobs
            jobs_completed += 1
            # Wake up all threads waiting on this condition.
            with job_cond:
                job_cond.notify_all()

        # --- Start of start function ---
        # File to show the status of the batch analyses
        statusfile = set_file_path(out+'.sta')
        with open(statusfile, mode='w', newline='\n') as file:
            file.write(f'{"Date":^10} {"Time":^8}  {"Batch":>6}  '
                       f'{"Iterations":>8}  {"Total":>6}  {"Empty":>14}  '
                       f'{"Incomplete":>14}  {"Completed":>14}\n')
            file.write(f'{now()}  {"0":>6}  {"0":>10}  {"0":>6}  {"0":>6}  '
                       f'{"0":>6}  {"0":>6}   {"0":>6}  {"0":>6} {"0":>6}')
        prv('statusfile initial size', statusfile.stat().st_size)
        # File to show the predictions of the batch analyses
        resultfile = set_file_path(out+'.csv')
        # Initialise class attribute
        self.debug = _init_debug()
        # Setup NumPy SeedSequence to produce repeatable pseudo-random numbers
        # across multiple processes.
        seedsequence = SeedSequence(_init_seed())
        # Setup batch counter
        batchcounter = count(start=0, step=1)
        batch = next(batchcounter)
        # Setup concurrent.futures variables
        running_futures = set()
        jobs_completed = 0
        percentage_completed = 0.0
        job_cond = threading.Condition()
        # This threading event flag is initially false.
        end_concurrent_execution = threading.Event()
        submit_job = threading.Event()
        submit_job.set()
        # Notes:
        # - The main cpu core "MainThread" will submit jobs to different cpus
        #   cores.
        # - The inner while-loop ensures continuous job submission.
        # - A "submit_job" event flag is used to ensure the number of jobs
        #   submitted do not exceed the max. number of cpus core that is to be
        #   used.
        # - When a job is completed, the main cpu core "Thread-1" will post-
        #   process the job's future.
        # - job_cond.wait() and job_cond.notify_all() are used to exit the job
        #   submission loop and prevent a deadlock scenario between the
        #   activities of the "MainThread" and "Thread-1".
        # - When the event flag end_concurrent_execution value is True, the
        #   concurrent submission of jobs to the cpu core pool will cease.
        #   executor will await all running_futures to be done before shutting
        #   down itself.
        # Start concurrent execution
        print(f'\n*** Start concurrent analyses of COVID19 confirmed case data. ***')
        with cf.ProcessPoolExecutor(max_workers=cpus) as executor:
            while True:  # job condition loop
                while percentage_completed < completion_target:  # job submission loop
                    # Submit job
                    job_fn = self.analyse
                    job_fn_args = _get_job_fn_args()
                    future = executor.submit(job_fn, *job_fn_args)
                    future.add_done_callback(_on_complete)
                    running_futures.add(future)
                    # Show submitted jobs in terminal
                    # stats = f'{now()}  cpus={cpus}  rf={len(running_futures)}  '\
                    #        f'batch={batch}  iterations={iterations*(batch)}'
                    # print(stats)
                    # Define event condition to pause job loop
                    if len(running_futures) >= cpus:
                        submit_job.clear()
                    # Allow job loop to proceed only when the event "submit_job"
                    # has a True value.
                    submit_job.wait()
                # Outside the job submission loop, here we check the Wait to be notified
                try:
                    with job_cond:
                        job_cond.wait()
                    # releases the underlying lock and then blocks until it is
                    # awakened by notify_all().
                    # call for the same condition variable in another thread.
                except KeyboardInterrupt:
                    print(f'\n## DETECTED {len(running_futures)} INCOMPLETED FUTURES ##')
                    print(f' ...ATTEMPTING TO CANCEL THESE FUTURES')
                    print(f' ...AFTER THESE FUTURES COMPLETES, PRESS CTRL+C MULTIPLE TIMES UNTIL SCRIPT ENDS.\n')
                    for future in running_futures:
                        _ = future.cancel()
                    # Ensure concurrent.futures.executor jobs really do finish.
                    _ = cf.wait(running_futures, timeout=None)
                # Condition to terminate the continuous streaming of job
                # submissions.
                if end_concurrent_execution.is_set():
                    print('\n*** Completed concurrent analyses of COVID19 '
                          'confirmed case data.  ***')
                    break
                if percentage_completed <= completion_target:
                    for future in running_futures:
                        _ = future.cancel()
                    # Ensure concurrent.futures.executor jobs really do finish.
                    _ = cf.wait(running_futures, timeout=None)
                    print(f'\n*** Satisfied the {completion_target}% '
                          f'Completion Target. ***')
                    break
        print(f'    Results are written in:\n {resultfile}.')
        # print( '\n*** Concurrent Analyses Ended ***' )

    def analyse(self, iterations, childseed, batch, days, ncases, mu_scenarios,
                sigma_scenarios, covid_earliest_day, factor, quota):
        '''Function to analyse COVID19 confirmed case data.

        Arguments:
         iterations - an integer; number of iterations to be performed
         batch - an integer; batch number
         rng - a numpy.random.default_rng()
         days - a list of day numbers
         ncases - number of daily cases
         mu_scenarios - a list of average-daily-confirmation-period to be
                        considered
         sigma_scenarios - a list of standard deviations for each mu_scenarios
         covid_earliest_day - the day number of the 1st imported case(s).

        A 3D NumPy array is used to store the no. of iterations, number of mu
        scenarios, number of daily cases

        '''
        def _get_norm_rsv(mus, sigmas, size):
            nonlocal childseed
            '''Returns a numpy.ndarray with the random variate of a normal
            distribution.

            Here, it is used to obtain a normal-distribution of the average
            number of day(s) taken to detect COVID19 after SARS-CoV2 infection
            and their standard deviations.

            A descendant of childseed is spawned each time this function is
            called. This descendant seed is used to generate an independent
            (with very high probability) streams of random number Generator
            that uses NumPy's PCG64DXSM BitGenerator.
            '''
            grandseed = childseed[0].spawn(1)
            #rng = default_rng(grandseed[0])
            rng = Generator(PCG64DXSM(grandseed[0]))
            cast_days = np.int16(
                norm.rvs(loc=mus, scale=sigmas, size=size, random_state=rng))
            return cast_days

        def _get_histogram(a, bins=10, range=None, normed=None, weights=None,
                           density=None):
            '''Returns the values of a histogram.'''
            hist, _ = np.histogram(a, bins=bins, range=range, normed=normed,
                                   weights=weights, density=density)
            return hist

        def _fix_prediction_ends(debug, label, predictions, origin, samples_mu,
                                 samples_sigma):
            '''Force predictions to be within bounds (1st imported COVID19 case
            day and last day of Local COVID-19 data.'''
            exceeds = np.where((predictions < covid_earliest_day) | (predictions > locallastday))
            if debug:
                predictions_mins = np.min(predictions, axis=2)
                predictions_maxs = np.max(predictions, axis=2)
                npri('predictions_mins', predictions_mins)
                npri('predictions_maxs', predictions_maxs)
                pri('exceeds', exceeds)
            if exceeds[0].shape[0] != 0:
                for r, c, d in zip(*exceeds):
                    while True:
                        #print( r, c, samples_mu[r,c,d], samples_sigma[r,c,d], origin[r,c,d] )
                        cast_day = _get_norm_rsv(samples_mu[r, c, d], samples_sigma[r, c, d], 1)
                        if label in 'backcast':
                            day = origin[r, c, d] - cast_day
                        elif label in 'forecast':
                            day = origin[r, c, d] - cast_day
                        #print( 'cast_day, day, covid_earliest_day', *cast_day, *day, covid_earliest_day )
                        if day >= covid_earliest_day and day <= locallastday:
                            predictions[r, c, d] = day[0]
                            break
            if debug:
                predictions_mins = np.min(predictions, axis=2)
                predictions_maxs = np.max(predictions, axis=2)
                npri('predictions_mins', predictions_mins)
                npri('predictions_maxs', predictions_maxs)
            return predictions

        def _select_mus(covid_days_bins, covid_days_ncases, ncases,
                        mu_scenarios, factor):
            # Make the size of ncases to be the same as size of
            # covid_days_ncases axis2 by padding a list of zeros at the start
            # of ncases.
            covid0days_index = np.where(covid_days_bins == 0)[0][0]
            lncases = np.insert(ncases,
                                np.zeros(covid0days_index, dtype=np.int16), 0)
            # print( 'covid0days_index', covid0days_index )
            # npri( 'lncases', lncases )
            results = {}
            for i in range(covid_days_bins[:-1].size):
                empirical = lncases[i]
                predicted = covid_days_ncases[:, i]
                # if empirical in predicted.tolist() :
                if (np.where(predicted >= factor*empirical)[0].any() &
                        np.where(predicted <= empirical)[0].any()):
                    same = np.where(predicted == empirical)[0]
                    selected_mus = [mu_scenarios[n] for n in same]
                    results.update({covid_days_bins[i]: selected_mus})
                    # print( f'{i} empirical={empirical}  predicted={predicted} '\
                    #       f'same={same} selected_mus={selected_mus}' )
            # pri( 'results', results )
            return results

        # --- start of analyse function ----
        # Show job batch information in terminal
        stats = f'{now()}  {threading.current_thread().getName()}  '\
                f'PID={os.getpid()}  batch={batch}  '\
                f'iterations={iterations*(batch)}'
        print(stats)
        # local variable
        locallastday = days[-1]
        # Convert array on daily mu and sigmas and days to every COVID19 case
        axis0 = iterations
        axis1 = mu_scenarios.size
        axis2 = np.sum(ncases).tolist()
        pshape = (axis0, axis1, axis2)
        samples_mu = np.broadcast_to(mu_scenarios[None, :, None], pshape)
        samples_sigma = np.broadcast_to(sigma_scenarios[None, :, None], pshape)
        all_days = np.repeat(days, ncases)
        samples_days = np.broadcast_to(all_days[None, None, :], pshape)
        if self.debug:
            prv('axis0', axis0)
            prv('axis1', axis1)
            prv('axis2', axis2)
            pri('pshape', pshape)
            npri('samples_mu', samples_mu)
            npri('samples_sigma', samples_sigma)
            npri('samples_days', samples_days)
        # Backcast: get random predictions of when SARS-Cov2 infections
        #           occurred for the differenct mus scenarios for each COVID19
        #           case.
        backcast_days = _get_norm_rsv(samples_mu, samples_sigma, pshape)
        inf_days = samples_days - backcast_days
        if self.debug:
            npri('backcast_days', backcast_days, 2)
            npri('inf_days', inf_days, 2)
        inf_days = _fix_prediction_ends(self.debug, 'backcast', inf_days,
                                        samples_days, samples_mu,
                                        samples_sigma)
        if self.debug:
            npri('inf_days', inf_days, 2)
        # Get ncases of the predicted SARS-Cov2 infections
        # FOR DEBUGGING
        #inf_days_bins = np.arange( covid_earliest_day, locallastday+2, 1, dtype=np.int16 )
        #inf_days_ncases = np.apply_along_axis( _get_histogram, 1, inf_days, bins=inf_days_bins )
        #itotaldays = np.sum( inf_days_ncases, 1)
        # if self.atype in 'debug':
        #    npri( 'inf_days_bins', inf_days_bins )
        #    npri( 'inf_days_ncases', inf_days_ncases, 1 )
        #    npri( 'itotaldays', itotaldays )
        # Forecast: get random predictions of when COVID19 was confirmed
        #           for the differenct mus scenarios for each SARS-CoV2 case
        forecast_days = _get_norm_rsv(samples_mu, samples_sigma, pshape)
        covid_days = inf_days + forecast_days
        if self.debug:
            npri('forecast_days', forecast_days, 2)
            npri('covid_days', covid_days, 2)
        covid_days = _fix_prediction_ends(self.debug, 'forecast', covid_days,
                                          inf_days, samples_mu, samples_sigma)
        # Get ncases of the predicted COVID19 confirmed cases
        covid_days_bins = np.arange(covid_earliest_day, locallastday+2, 1,
                                    dtype=np.int16)
        covid_days_ncases = np.apply_along_axis(_get_histogram, 2, covid_days,
                                                bins=covid_days_bins)
        ctotaldays = np.sum(covid_days_ncases, 2)
        if self.debug:
            npri('covid_days', covid_days, 2)
            npri('covid_days_bins', covid_days_bins)
            npri('covid_days_ncases', covid_days_ncases, 2)
            npri('ctotaldays', ctotaldays)
        # Consolidate predicted mus for current batch
        batch_selected_mus = {i: [] for i in
                              range(covid_earliest_day, locallastday+1)}
        # a dict; {day number: a list of predicted mus}
        for n, each_iteration_covid_days_ncases in enumerate(covid_days_ncases):
            # For each day, check if the predicted number of COVID19 cases is
            # similar to empirical data. If so, select the corresponding mu for
            # that day.
            selected_mus = _select_mus(covid_days_bins,
                                       each_iteration_covid_days_ncases,
                                       ncases, mu_scenarios, factor)
            if self.debug:
                pri(f'{n} selected_mus', selected_mus)
            # Collate the selected mus for this batch if there isn
            for k, v in selected_mus.items():
                batch_selected_mus[k].extend(v)
        return (batch_selected_mus, batch, iterations)


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


def write_predicted_mus(days, dates, mus, headers, target='out.csv'):
    '''Function to write predicted mus into a .csv file.'''
    with open(target, mode='w', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile, delimiter=',', quoting=csv.QUOTE_NONE)
        # 1st row: write column labels
        writer.writerow(headers)
        # Subesquent rows: write values of days, dates, mus
        for i, j, k in zip(days, dates, mus,):
            if isinstance(k, list):
                writer.writerow([i, j, *k])
            else:
                writer.writerow([i, j, k])


def write_status(status, target='out.sta'):
    '''Function to write the progress of the prediction into file.'''
    with open(target, mode='a', newline='\n') as file:
        file.write(f'\n{status}')


def set_file_path(file):
    '''Return the full file name of file as a pathlib.Path object.'''
    source = Path.cwd()/Path(file)
    return source


def convert_str_to_datetime(date, dateformat="%d/%m/%Y", tzinfo=timezone.utc):
    refdate = datetime.strptime(date, dateformat)
    refdate = refdate.replace(tzinfo=tzinfo)
    return refdate


def now():
    '''Function returns currnet date and time in "dd/mm/yyyy hh:mm:ss" format.
    '''
    return datetime.now().strftime("%d/%m/%Y %H:%M:%S")


def initialize(covid19_csvfile, mu_min, mu_max, sigma_a, sigma_b,
               earliest_arrival_date_of_imported_cases=None):
    '''Function to initialize the arguments and class object needed by this
    module to start the prediction.'''

    def _get_covid19_case_data():
        '''Function to get Daily Imported and Local Confirmed COVID19 case data.

        Argument:
         covid19_csvfile - a ".csv" file; stores Daily Imported and Local
                                          Confirmed COVID19 case data
        Return:
         first_imported_case_date - a datetime.datetime class object;
                                    format "2020-01-23 00:00:00+00:00"
         local - a Local_COVID19_Parameters class object;
         '''
        # Read-in Daily Imported and Local Confirmed COVID19 cases and store
        #  them as unique dictionaries.
        imported_covid_dict, local_covid_dict = _get_daily_cases(datafile)
        # Get earliest possible date for local COVID19.
        first_imported_case_date = min(imported_covid_dict.keys())
        if earliest_arrival_date_of_imported_cases is None:
            earliest_possible_local_covid_date = first_imported_case_date
        else:
            earliest_possible_local_covid_date = earliest_arrival_date_of_imported_cases
        prv('earliest_possible_local_covid_date', earliest_possible_local_covid_date)
        # Extract daily data from local_covid
        local = Local_COVID19_Parameters(local_covid_dict)
        return (local, earliest_possible_local_covid_date)

    def _get_daily_cases(csvfile):
        '''Extract the number of Daily Imported and Local Confirmed COVID19
           cases from "csvfile". Return corresponding dictionaries.'''
        imported = {}
        local = {}
        with open(csvfile, mode='r', encoding='utf-8') as csvfile:
            reader = csv.reader(csvfile)
            next(reader)  # Skip header row
            for row in reader:
                if row:
                    # print(f'row -> {type(row)} {row}')
                    date = convert_str_to_datetime(row[0],
                                                   dateformat='%Y-%m-%d')
                    imported[date] = int(row[1])
                    local[date] = int(row[2])
                    # print(date, imported, local)
        index = count()
        for k, v in imported.items():
            print(f'{next(index)} {k}  {v}  {local[k]}')
        return imported, local

    def _get_mu_sigma(mu_min, mu_max, a, b):
        '''Function to define the average number of days taken to confirm
        local covid19 cases and their corresponding standard deviations.

        Arguments:
         mu_min and mu_max - min. and max. average number of days taken to
                             locally confirm COVID19 after SARS-CoV2 infection.
         a,b - curve fitting constants for a quadratic mu-sigma relationship.

        Return:
         mu_sigma - a dictionary of the mu-sigma relationship
                    key = a range of the average number of days taken to
                          locally confirm COVID19 after SARS-CoV2 infection
                    value = standard deviations of the key
        '''
        mu = [i for i in range(mu_min, mu_max, 1)]
        mu_sigma = {x: a*x**2 + b*x for x in mu}
        pri('mu_sigma', mu_sigma)
        return mu_sigma

    # ---- start for initialize function ---
    earliest_arrival_date_of_imported_cases = convert_str_to_datetime(
        earliest_arrival_date_of_imported_cases,
        dateformat='%Y-%m-%d',
        tzinfo=timezone.utc)
    # Get covid19 confirmed case data.
    datafile = set_file_path(covid19_csvfile)
    local, earliest_covid_date = _get_covid19_case_data()
    # Define the average and standard deviation trend to be used to analyse
    # covid19 confirmed local cases.
    mu_sigma = _get_mu_sigma(mu_min, mu_max+1, sigma_a, sigma_b)
    # Initialize class object "Predict_COVID19_Confirmation_Periods"
    return Predict_COVID19_Confirmation_Periods(mu_sigma, local,
                                                earliest_covid_date)


def main():
    # Inputs needed for prediction
    covid19_csvfile = 'COVID19_epidemic_trends.csv'         # To be defined
    earliest_arrival_date_of_imported_cases = '2020-01-18'  # To be defined
    mu_min = 1            # To be defined
    mu_max = 18           # To be defined
    sigma_a = -0.008665   # To be defined
    sigma_b = 0.483888    # To be defined
    seed = 91023483923    # Random Seed1
    # seed = 66576402094  # Random Seed2
    # seed = 343090589475 # Random Seed3
    iterations = 150      # Need to be tuned to optimise computation
    cpus = 28             # Need to be tuned to optimise computation
    completion_target = 100.0
    debug = False         # Set to True for debug purposes
    # Define covid19_csvfile path
    repo_dir = Path(__file__).resolve().parents[2]
    covid_data = repo_dir / '4_Empirical_Data' / covid19_csvfile
    # Initialize Data and variables
    predict = initialize(covid_data, mu_min, mu_max, sigma_a, sigma_b,
                         earliest_arrival_date_of_imported_cases)
    # Start esimation
    predict.start(seed=seed, cpus=cpus, iterations=iterations, quota=300,
                  factor=1.0,  completion_target=completion_target,
                  debug=debug)
    print(f'#### COMPLETED completion_target of {completion_target}% ####')


if __name__ == "__main__":
    main()
