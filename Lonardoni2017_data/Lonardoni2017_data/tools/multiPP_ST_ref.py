import numpy as np
import os
import time
import psutil


###########################################################################################
###

def runFunctionsInParallel(listOf_FuncAndArgLists, kwargs=None, names=None, offsetsSeconds=None,
                           expectNonzeroExit=False, maxAtOnce=None, SZORDAN=1, sleepTime=.5):
    import time
    start = time.time()
    if not listOf_FuncAndArgLists:
        return []  # list of functions to run was empty.
    if SZORDAN:
        cpu_max = psutil.cpu_count() - 1
    else:
        cpu_max = maxAtOnce
    if offsetsSeconds is None:
        offsetsSeconds = 0

    # If no argument list is given, make one:
    N = len(listOf_FuncAndArgLists)
    listOf_FuncAndArgLists = [faal if isinstance(faal, list) else [faal, [], {}] for faal in listOf_FuncAndArgLists]
    listOf_FuncAndArgLists = [faal + [{}] if len(faal) == 2 else faal for faal in listOf_FuncAndArgLists]
    listOf_FuncAndArgLists = [faal + [[], {}] if len(faal) == 1 else faal for faal in listOf_FuncAndArgLists]



    from multiprocessing import Process, Queue

    if names is None:
        names = [None for _ in listOf_FuncAndArgLists]
    assert len(names) == len(listOf_FuncAndArgLists)

    def functionWrapper(fff, que, theArgs=None, kwargs=None,
                        delay=None):  # add a argument to function for assigning a queue
        os.nice(10)  # Add 10 to the niceness of this process (POSIX only)
        if delay:
            from time import sleep
            sleep(delay)
        if not kwargs:
            kwargs = {}
        if not theArgs:
            theArgs = []

        # print 'MULTIPROCESSING: Launching %s in parallel '%funcName
        returnVal = que.put(fff(*theArgs, **kwargs))
        # print 'MULTIPROCESSING: Finished %s in parallel! '%funcName
        return returnVal

    def reportStatus(sjobs, showmax, showsuccessful=10 ** 10, print_info=0):  # jobs):
        djobs = sjobs if showmax >= len(sjobs) else sjobs[:showmax]

        tableFormatString = '%s\t%' + str(max([len(job.name) for job in djobs])) + 's:\t%9s\t%s\t%s\t%s'
        if print_info: print('\n' + '-' * 75 + '\n' + tableFormatString % ('alive?', 'Job', 'exit code', 'Full', 'Empty', 'Func',) + '\n' + '-' * 75)
        # Check that we aren't going to show more *successfully finished* jobs than we're allowed: Find index of nth-last successful one. That is, if the limit binds, we should show the latest N=showsuccessful ones only.
        isucc = [ii for ii, job in enumerate(djobs) if job.exitcode == 0]
        earliestSuccess = -1 if len(isucc) < showsuccessful else isucc[::-1][showsuccessful - 1]
        if print_info: print('\n'.join([tableFormatString % (job.is_alive() * 'Yes:', job.name, job.exitcode, queues[iii].full(), queues[iii].empty(),'(built-in function)' if not hasattr(listOf_FuncAndArgLists[ii][0], 'func_name') else listOf_FuncAndArgLists[ii][0].func_name) for ii, job in enumerate(djobs) if job.exitcode != 0 or ii >= earliestSuccess]))

        if len(isucc) > showsuccessful:
            if print_info: print('%d other jobs finished successfully.' % (len(isucc) - showsuccessful))
        if len(sjobs) > len(djobs):
            if print_info: print('%d more jobs waiting for their turn to start...' % (len(sjobs) - len(djobs)))
        if print_info: print('-' * 75 + '\n')
        return [job.exitcode for ii, job in enumerate(sjobs)]

    def emptyQueues():  # jobs,queues,gotQueues):
        for ii, job in enumerate(jobs):
            if not queues[ii].empty():
                if ii in gotQueues:
                    gotQueues[ii] += queues[ii].get()
                else:
                    gotQueues[ii] = queues[ii].get()

    def fillQueues(check):  # jobs,queues,gotQueues):
        flag = (int(10 * (1 - psutil.cpu_percent(check) * .01)) * .8)
        for ii, job in enumerate(jobs):
            #            print ii
            if not queues[ii].empty():

                gotQueues.append(queues[ii].get())
                if istart > cpu_max:
                    jobs[ii].join()

                if any(listOf_FuncAndArgLists) and ii <= maxAtOnce:
                    funcArgs = listOf_FuncAndArgLists.pop()


                    jobs[ii] = Process(target=functionWrapper, args=[funcArgs[0], queues[ii]],
                                       kwargs={'theArgs': funcArgs[1], 'kwargs': funcArgs[2], 'delay': delays[ii]},
                                       name=names[ii])
                    jobs[ii].start()
            elif SZORDAN and not job.is_alive() and ii > maxAtOnce and flag > 0:
                if any(listOf_FuncAndArgLists):
                    funcArgs = listOf_FuncAndArgLists.pop()
                    jobs[ii] = Process(target=functionWrapper, args=[funcArgs[0], queues[ii]],
                                       kwargs={'theArgs': funcArgs[1], 'kwargs': funcArgs[2], 'delay': delays[ii]},
                                       name=names[ii])
                    jobs[ii].start()
                    flag += -1

    queues = [Queue() for _ in xrange(maxAtOnce)]  # create a queue object for each function
    delays = list(np.arange(cpu_max) * offsetsSeconds)

    jobs = [Process(target=functionWrapper, args=[funcArgs[0], queues[iii]],
                    kwargs={'theArgs': funcArgs[1], 'kwargs': funcArgs[2], 'delay': delays[iii]}, name=names[iii]) for
            iii, funcArgs in enumerate(listOf_FuncAndArgLists[:maxAtOnce])]
    listOf_FuncAndArgLists = listOf_FuncAndArgLists[maxAtOnce:]
    #    print len(listOf_FuncAndArgLists)
    cpu_max = min(cpu_max, len(listOf_FuncAndArgLists))
    istart = cpu_max if cpu_max < len(jobs) else maxAtOnce
    istart = 0
    while len(jobs) < cpu_max:
        jobs.append(Process(target=fake_process))
        queues.append(Queue())
    for job in jobs: job.start()  # Launch them all
    #    print len(jobs),cpu_max

    import time

    timeElapsed = 0
    gotQueues = []


    """ Now, wait for all the jobs to finish.  Allow for everything to finish quickly, at the beginning.
    """
    while istart < N or any([jj.is_alive() for jj in jobs]):
        #        if istart==4096:
        #            print [jj.is_alive() for jj in jobs]

        #        if timeElapsed>0:
        #            time.sleep(sleepTime) # Wait a while before next update. Slow down updates for really long runs.
        timeElapsed += sleepTime
        # Add any extra jobs needed to reach the maximum allowed:
        #        while istart<len(jobs) and sum([jj.is_alive() for jj in jobs]) < maxAtOnce:
        #            jobs[istart].start()
        #            istart+=1
        #            timeElapse=.01
        #        reportStatus(jobs,istart,showFinished) #istart)#jobs)
        fillQueues(sleepTime)
        istart = len(gotQueues)
    # print istart
    #        emptyQueues()#jobs,queues,gotQueues)
    #        sys.stdout.write("\r %d more jobs waiting for their turn to start..." % (N-istart))

    #        sys.stdout.flush()

    for job in jobs: job.join()  # Wait for them all to finish... Hm, Is this needed to get at the Queues?

    # And now, collect any remaining buffered outputs (queues):
    emptyQueues()
    #    print jobs
    for ii, job in enumerate(jobs):
        #        print ii,
        if len(gotQueues[ii]) == 0:
            gotQueues[ii] = None
    print " "
    # Give final report of exit statuses?
    success = reportStatus(jobs, np.inf)
    # names=[names[iii] if names[iii] is not None else fff[0].func_name for iii,fff in enumerate(listOf_FuncAndArgLists)]

    if any(success):
        print('MULTIPROCESSING: Parallel processing batch set did not ALL succeed successfully ')
        assert expectNonzeroExit  # one of the functions you called failed.
        return False
    else:
        print('MULTIPROCESSING: Apparent success of all functions')
        print ('MULTIPROCESSING: %d s elapsed for computation' % (time.time() - start))
    return success, [gotQueues[ii] for ii in range(len(gotQueues))]


def fake_process(a=1, b=2, c=3, d=4, e=5, f=6, g=7, h=8, i=9):
    time.sleep(1)
    return -123.45678