#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Extract data for plotting
"""

import os
import sys
import numpy as np

import glob

#import pylab as pl
#import matplotlib as mpl
#from matplotlib.backends.backend_pdf import PdfPages

#from blist import sortedset

def die_with_usage():
    """ HELP MENU """
    print 'main.py'
    print ''
    print 'USAGE'
    print '   ./main.py <input>'
    print 'PARAMS'
    print '          output - input file name'
    sys.exit(0)

if __name__ == '__main__':

    ## help menu
    if len(sys.argv) < 1:
        die_with_usage()

    ## params
    pathName = sys.argv[1]
    outFName = sys.argv[2]

    if os.path.isfile(outFName):
        print 'WARNING:',outFName,'already exists.'
        print 'Overwriting file................'
        #sys.exit(0)
    
    fout = open(outFName, 'w')
    #fout.write("#       core    g_x_g_y_g_z_v#      BPF_BPB     FC      AF      AC      IDX_xXyXz    Time(m)        Thuput(m)        Create(m)      rst(m)   brst(m)  hz(m)   agg(m)  io(m)   Time(s) Thuput(s)       Create(s)       rst(st) brst(st) hz(st) agg(st)  io(st)\n")
    
    fout.write("#       core    g_x_g_y_g_z_v#  Bx*By*Bz      BPF_BPB     FC      AF      AC      x*y*z    Time(m)        Thpt(m)       Thpt(max)        Create(m)      rst(m)   brst(m)  hz(m)   agg(m)  io(m)\n")
    
    for inpFName in sorted(glob.glob(os.path.join(pathName,'*'))):
    
        print inpFName
        finp = open(inpFName, 'r')

        ## start loop
        skip_first_1 = 0
        skip_first_2 = 0
        skip_first_3 = 0

        
        time = []
        throughput = []
        file_create_time = []
        rst_time = []
        b_rst_time = []
        hz_time = []
        agg_time = []
        io_time = []
        brst_x = 1;
        brst_y = 1;
        brst_z = 1;
        
        iters = 0
        for line in finp:
            if 'Blocks Per File' in line:
                words = line.split()
                bpf = words[3]
                bpb = words[7]
                fc = words[10]
                af = words[13]
                ac = words[16]
                
            if 'Blocks Restructuring Size' in line:
                words = line.split()
                brst_x = words[3]
                brst_y = words[4]
                brst_z = words[5]
                                
            if 'Global Data' in line:
                words = line.split()
                cores = words[1]
                gl_x = words[4]
                gl_y = words[5]
                gl_z = words[6]
                variables = words[8]
                idx_count = words[11]
                idx_count_x = words[13]
                idx_count_y = words[15]
                idx_count_z = words[17]
                
            if 'File Create time' in line:
                if skip_first_1 == 0:
                    skip_first_1 = 1
                else:
                    words = line.split()
                    file_create_time.append(float(words[6]))
                    
            if 'Write time' in line:
                if skip_first_2 == 0:
                    skip_first_2 = 1
                else:
                    words = line.split()
                    rst_time.append(float(words[11]))
                    b_rst_time.append(float(words[13]))
                    hz_time.append(float(words[15]))
                    agg_time.append(float(words[17]))
                    io_time.append(float(words[19]))
                    
            
            if 'Time Taken' in line:
                if skip_first_3 == 0:
                    skip_first_3 = 1
                else:
                    words = line.split()
                    time.append(float(words[2]))
                    throughput.append(float(words[5]))                    
                    #print time
                    #print throughput
                    iters = iters + 1

        meanTime = np.asarray(time).mean()
        stddevTime = np.asarray(time).std()
        maxTime = np.asarray(time).max()
        
        meanfileCreate = np.asarray(file_create_time).mean()
        stddevfileCreate = np.asarray(file_create_time).std()
        maxfileCreate = np.asarray(file_create_time).max()
        
        meanRST = np.asarray(rst_time).mean()
        stddevRST = np.asarray(rst_time).std()
        maxRST = np.asarray(rst_time).max()
        
        meanBRST = np.asarray(b_rst_time).mean()
        stddevBRST = np.asarray(b_rst_time).std()
        maxBRST = np.asarray(b_rst_time).max()
        
        meanHZ = np.asarray(hz_time).mean()
        stddevHZ = np.asarray(hz_time).std()
        maxHZ = np.asarray(hz_time).max()
        
        meanAGG = np.asarray(agg_time).mean()
        stddevAGG = np.asarray(agg_time).std()
        maxAGG = np.asarray(agg_time).max()
        
        meanIO = np.asarray(io_time).mean()
        stddevIO = np.asarray(io_time).std()
        maxIO = np.asarray(io_time).max()
        
        meanThroughput = np.asarray(throughput).mean()
        stddevThroughput = np.asarray(throughput).std()
        maxThroughput = np.asarray(throughput).max()
        
        #fout.write('%d  %s  %sx%sx%sx%s   %sx%s      %s      %s      %s %sx%sx%s        %f      %f      %f      %f      %f      %f      %f      %f      %f      %f      %f      %f      %f      %f      %f      %f \n' % (iters, cores, gl_x, gl_y, gl_z, variables, bpf, bpb, fc, af, ac, idx_count_x, idx_count_y, idx_count_z, meanTime, meanThroughput, meanfileCreate, meanRST, meanBRST, meanHZ, meanAGG , meanIO, stddevTime, stddevThroughput, stddevfileCreate, stddevRST, stddevBRST, stddevHZ, stddevAGG, stddevIO))
        
        fout.write('%d  %s  %sx%sx%sx%s %sx%sx%s   %sx%s      %s      %s      %s %sx%sx%s        %f      %f      %f      %f      %f      %f      %f      %f      %f\n' % (iters, cores, gl_x, gl_y, gl_z, variables, brst_x, brst_y, brst_z, bpf, bpb, fc, af, ac, idx_count_x, idx_count_y, idx_count_z, meanTime, meanThroughput, maxThroughput, meanfileCreate, meanRST, meanBRST, meanHZ, meanAGG , meanIO))

        finp.close()
    
    fout.close()