# -*- coding: utf-8 -*-

# This file is part of OpenDrift.
#
# OpenDrift is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2
#
# OpenDrift is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with OpenDrift.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2017, André R. Brodtkorb, SINTEF Digital

import multiprocessing
import Queue
from timeit import default_timer as timer

import glob
import os
import sys


def run_test(out_q, path):
    import imp
    
    ##########################################
    from opendrift.models.basemodel import OpenDriftSimulation
    from opendrift.models.opendrift3D import OpenDrift3DSimulation
    from opendrift.models.openoil3D import OpenOil3D
    
    # Monkey patch opendrift to avoid plotting
    def noplot(*args, **kwargs):
        pass
    OpenDriftSimulation.plot = noplot
    OpenDriftSimulation.animation = noplot
    OpenDriftSimulation.plot_property = noplot
    OpenDriftSimulation.animation_profile = noplot
    OpenDrift3DSimulation.plot_vertical_distribution = noplot
    OpenOil3D.plot_oil_budget = noplot
    
    perflog = []
    
    # Monkey patch opendrift to save performance log
    orig_perf = OpenDriftSimulation.performance
    def saveperf(self):
        log = orig_perf(self)
        perflog.append(log)
        return log
    OpenDriftSimulation.performance = saveperf
    
    ##########################################
    
    # Override stdout
    filename = os.path.basename(path)
    sys.stdout = sys.stderr = open(str(filename) + ".out", "w")
    
    # Import module - execute script
    start = timer()
    try:
        temp_module = imp.load_source('module.name', path)
    except:
        print("FAILED")
    end = timer()
    
    # Write out total time
    duration = (end-start)
    out_q.put([ path, duration, perflog] )
    
    
def main():
    out_q = multiprocessing.Queue(1)
    
    results = []
    
    examples_path = os.path.dirname(os.path.realpath(__file__))
    files = glob.glob(os.path.join(examples_path, 'example_*.py'))
    

    # Set backend of matplotlib 
    #import matplotlib
    #import matplotlib.pyplot as plt
    #plt.switch_backend('agg')
    
    
    def printResults(result):
        print('********************************************')
        print(result[0])
        print('Total time: ' + str(result[1]))
        print('********************************************')
        for log in result[2]:
            print('==============================================')
            print(log)
            print('==============================================')
        print('********************************************')
    
    start = timer()
    for idx, file in enumerate(files):
            
        print('FILE ' + str(idx) + '/' + str(len(files)) + ': ' + str(file))
        proc = multiprocessing.Process(target=run_test, args=(out_q, file))
        print('Starting ' + str(proc))
        proc.start()
        
        try:
            var = out_q.get(timeout=600)
            results.append(var)
        except Queue.Empty:
            print("Failed - No result after timeout!")
            exit
        #printResults(var)
        
        proc.join(timeout=10)
        print('Joined ' + str(proc))
        
        print('Elapsed: ' + str(timer() - start))
        
        
    for result in results:
        printResults(result)
        
if __name__ == '__main__':
    main()
