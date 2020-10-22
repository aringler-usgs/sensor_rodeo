#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

day = 289
flag = 'Mean Ratio'


if flag == 'Orient':
    f = open('Results_'+ str(day) + '.csv','r')
    vals12, vals23, vals13, times = [], [], [], []
    for line in f:
        if flag in line:
            line = line.split(', ')
            if '1 2' in line[0]:
                vals12.append(float(line[3]))
                times.append(float(line[1]) + float(line[2])/24.)
            if '2 3' in line[0]:
                vals23.append(float(line[3]))
            if '1 3' in line[0]:
                vals13.append(float(line[3]))

    f.close()

    fig = plt.figure(1,figsize=(8,8))
    plt.subplot(3,1,1)
    plt.plot(times,vals12,'.', label='1 relative to 2')
    plt.legend(loc=2)
    plt.ylim((-5,5))
    plt.subplot(3,1,2)
    plt.plot(times,vals23,'.', label='2 relative to 3')
    plt.legend(loc=2)
    plt.ylim((-5,5))
    plt.subplot(3,1,3)
    plt.plot(times,vals13,'.', label='1 relative to 3')
    plt.legend(loc=2)
    plt.ylim((-5,5))
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




