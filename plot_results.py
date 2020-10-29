#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times')
mpl.rc('text', usetex=True)
mpl.rc('font',size=16)


locs = ['West Wault', 'Cross-Tunnel', 'East-Tunnel']
days = [282, 289, 296]
idx = 2

flag = 'Self-Noise'
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
    print(vals12)
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


elif flag == 'Mean Ratio':
    fig = plt.figure(1,figsize=(8,8))
    for idx, comp in enumerate(['Z', '1', '2']):
        f = open('Results_'+ str(day) + '.csv','r')
        vals12, vals23, vals13, times = [], [], [], []
        std12, std23, std13 = [], [], []
        for line in f:
            if flag in line:
                line = line.split(', ')
                if '1 2' in line[0]:
                    vals12.append(float(line[4]))
                    std12.append(float(line[5]))
                    times.append(float(line[2]) + float(line[3])/24.)
                if '2 3' in line[0]:
                    vals23.append(float(line[4]))
                    std23.append(float(line[5]))
                if '1 3' in line[0]:
                    vals13.append(float(line[4]))
                    std13.append(float(line[5]))
        f.close()
        plt.subplot(3,1, idx + 1)
        plt.errorbar(times, vals12, yerr=std12, marker='.', linestyle='',  label='STS-6 relative to T-360')
        plt.errorbar(times, vals23, yerr= std23, marker='.', linestyle='', label='T-360 relative to STS-2.5')
        plt.errorbar(times, vals13,yerr=std13, marker='.', linestyle='',  label='STS-6 relative to STS-2.5')
        plt.text(min(times), 1.007, 'BH' + comp)
        plt.xlim((min(times)-.1, max(times)+.1))

        plt.legend(loc='lower center', ncol=3, fontsize=10)
        if idx == 1:
            plt.ylabel('Mean Ratio')
        if idx == 0:
            plt.title(loc)
        plt.ylim((0.99,1.01))
    plt.xlabel('Time (doy)')
    plt.savefig(flag.replace(' ','') + '_' + str(day) + '.png', format='PNG')

else:
    # Need to fix the overplotting
    fig = plt.figure(1,figsize=(8,8))
    for idx, comp in enumerate(['Z', '1', '2']):
        f = open('Results_'+ str(day) + '.csv','r')
        vals1, vals2, vals3, times = [], [], [], []
        std1, std2, std3 = [], [], []
        for line in f:
            if flag in line:
                line = line.split(', ')

                if ('1' in line[1]) and (comp in line[2]):
                    vals1.append(float(line[5]))
                    std1.append(float(line[6]))
                    times.append(float(line[3]) + float(line[4])/24.)
                if ('2' in line[1]) and (comp in line[2]):
                    vals2.append(float(line[5]))
                    std2.append(float(line[6]))
                if ('3' in line[1]) and (comp in line[2]):
                    vals3.append(float(line[5]))
                    std3.append(float(line[6]))
        f.close()
        plt.subplot(3,1, idx + 1)
        plt.errorbar(times, vals1, yerr=std1, marker='.', linestyle='',  label='STS-6')
        plt.errorbar(times, vals2, yerr= std2, marker='.', linestyle='', label='T-360')
        plt.errorbar(times, vals3,yerr=std3, marker='.', linestyle='',  label='STS-2.5')
        if flag == 'Power':
            plt.text(min(times), -145, 'BH' + comp)
        else:
            plt.text(min(times), -165, 'BH' + comp)
        plt.xlim((min(times)-.1, max(times)+.1))
        plt.legend(loc='lower center', ncol=3, fontsize=10)
        if idx == 1:
            plt.ylabel('dB relative to 1 $(m/s^2)^2/Hz$')
        if idx == 0:
            plt.title(loc)
        if flag == 'Power':
            plt.ylim((-167,-140))
        else:
            plt.ylim((-182,-160.))
    plt.xlabel('Time (doy)')
    plt.savefig(flag.replace(' ','') + '_' + str(day) + '.png', format='PNG')

