#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times')
mpl.rc('text', usetex=True)
mpl.rc('font',size=16)




locs = ['West Wault', 'Cross-Tunnel']
days = [282, 289]
idx = 0

flag = 'Orient'
day = days[idx]
loc = locs[idx]

if flag == 'Orient':
    f = open('Results_'+ str(day) + '.csv','r')
    vals12, vals23, vals13, times = [], [], [], []
    std12, std23, std13 = [], [], []
    for line in f:
        if flag in line:
            line = line.split(', ')
            val = float(line[3])
            if val > 6:
                val = 6
            if val < -6:
                val = -6
            if '1 2' in line[0]:
                vals12.append(val)
                times.append(float(line[1]) + float(line[2])/24.)
                std12.append(float(line[4]))
            if '2 3' in line[0]:
                vals23.append(val)
                std23.append(float(line[4]))
            if '1 3' in line[0]:
                vals13.append(val)
                std13.append(float(line[4]))
    f.close()
    print(std12)
    fig = plt.figure(1,figsize=(8,8))
    plt.subplot(3,1,1)
    plt.title(loc)
    plt.plot(times,vals12,marker='.', linestyle='', label='STS-6 relative to T-360')
    plt.legend(loc=4)
    plt.ylim((-6,6))
    plt.xlim((min(times), max(times)))
    plt.subplot(3,1,2)
    plt.plot(times,vals23, marker='.', linestyle='', label='T-360 relative to STS-2.5')
    plt.legend(loc=4)
    plt.ylabel('Relative Orientation (degree)')
    plt.ylim((-6,6))
    plt.xlim((min(times), max(times)))
    plt.subplot(3,1,3)
    plt.plot(times,vals13,marker='.', linestyle='', label='STS-6 relative to STS-2.5')
    plt.legend(loc=4)
    plt.ylim((-6,6))
    plt.xlim((min(times), max(times)))
    plt.savefig(flag + '_' + str(day) + '.png', format='PNG')

if flag == 'Mean Ratio':
    fig = plt.figure(1,figsize=(8,8))
    for idx, comp in enumerate(['Z', '1', '2']):
        f = open('Results_'+ str(day) + '.csv','r')
        vals12, vals23, vals13, times = [], [], [], []
        for line in f:
            if flag in line:
                line = line.split(', ')
                if '1 2' in line[0]:
                    vals12.append(float(line[4]))
                    times.append(float(line[2]) + float(line[3])/24.)
                if '2 3' in line[0]:
                    vals23.append(float(line[4]))
                if '1 3' in line[0]:
                    vals13.append(float(line[4]))
        f.close()
        plt.subplot(3,1, idx + 1)
        plt.plot(times, vals12, '.', label='1 relative to 2')
        plt.plot(times, vals23, '.', label='2 relative to 3')
        plt.plot(times, vals13, '.', label='1 relative to 3')
        plt.text(min(times), 1.005, 'BH' + comp)
        plt.legend()
        if idx == 1:
            plt.ylabel('Mean Ratio')
        plt.ylim((0.99,1.01))
plt.xlabel('Time (doy)')
plt.savefig(flag.replace(' ','') + '_' + str(day) + '.png', format='PNG')




