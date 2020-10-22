#!/usr/bin/env python
import glob
import numpy as np
import matplotlib.pyplot as plt
from obspy.core import UTCDateTime, read, Stream
from obspy.core.inventory import read_inventory
from matplotlib.mlab import csd
from scipy.optimize import root
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm



import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times')
mpl.rc('text', usetex=True)
mpl.rc('font',size=16)

debug, plot = True, True
stimes, etimes = [], []
length, overlap, windows = 2**14, 2**8, 4*60*60

# West vault firs times
stimes.append(UTCDateTime('2020-282T21:00:00'))
etimes.append(UTCDateTime('2020-289T16:00:00'))

# Xtunnel vault need to adjust end time
stimes.append(UTCDateTime('2020-289T17:45:00'))
etimes.append(UTCDateTime('2020-296T00:00:00'))

def cp(tr1,tr2):
    cpval,fre = csd(tr1.data, tr2.data, NFFT=length,
                    Fs=tr1.stats.sampling_rate,
                    noverlap=overlap, scale_by_freq=True)
    return cpval[1:], fre[1:]

def plot_psd_noise(psd1, n11, psd2, n22, psd3, n33, per, st_chan):
    fig = plt.figure(1, figsize=(9,9))
    plt.semilogx(per,psd1, label='PSD ' + (st_chan[0].id).replace('.', ' '))
    plt.semilogx(per,n11, label='Noise ' + (st_chan[2].id).replace('.', ' '))
    plt.semilogx(per,psd2, label='PSD ' + (st_chan[1].id).replace('.', ' '))
    plt.semilogx(per,n22, label='Noise ' + (st_chan[2].id).replace('.', ' '))
    plt.semilogx(per,psd3, label='PSD ' + (st_chan[2].id).replace('.', ' '))
    plt.semilogx(per,n33, label='Noise ' + (st_chan[2].id).replace('.', ' '))
    per2, nlnm = get_nlnm()
    per2, nhnm = get_nhnm()
    plt.semilogx(per2, nlnm, color='k', linewidth=2)
    plt.semilogx(per2, nhnm, color='k', linewidth=2, label='NLNM/NHNM')
    plt.xlabel('Period (s)')
    plt.ylabel('PSD (dB rel. 1 $(m/s^2)^2/Hz$)', fontsize=16)
    plt.legend(ncol=2)
    plt.xlim((1./20., 300.))
    plt.title('PSD '+ str(st_chan[0].stats.starttime.julday).zfill(3) + ' ' +
        str(st_chan[0].stats.starttime.hour).zfill(2) + ':' +
        str(st_chan[0].stats.starttime.minute).zfill(2) + ' ' + st_chan[0].stats.component)
    plt.savefig('PSD_' + str(st_chan[0].stats.starttime.julday).zfill(3) + '_' +
        str(st_chan[0].stats.starttime.hour).zfill(2) + '.png')

    plt.clf()


    return

def write_ratios(psd1, psd2, psd3, per, f, time, chan):
    mb1 = psd1[(4.0 <= per) & (per <= 8.0)]
    mb2 = psd2[(4.0 <= per) & (per <= 8.0)]
    mb3 = psd3[(4.0 <= per) & (per <= 8.0)]
    diff1 = np.average(mb1/mb2)
    diff2 = np.average(mb2/mb3)
    diff3 = np.average(mb1/mb3)
    std1 = np.std(mb1/mb2)
    std2 = np.std(mb2 /mb3)
    std3 = np.std(mb1 /mb3)
    f.write('Mean Ratio 1 2, ' + chan + ', ' + str(time.julday).zfill(3) + ', ' + 
            str(time.hour).zfill(2) + ', ' + str(diff1.real) + ', ' + str(std1) + '\n')
    f.write('Mean Ratio 1 3, ' + chan + ', ' + str(time.julday).zfill(3) + ', ' + 
            str(time.hour).zfill(2) + ', ' + str(diff3.real) + ', ' + str(std3) + '\n')
    f.write('Mean Ratio 2 3, ' + chan + ', ' + str(time.julday).zfill(3) + ', ' + 
            str(time.hour).zfill(2) + ', ' + str(diff2.real) + ', ' + str(std2) + '\n')
    return

def write_noise(psd1, psd2, psd3, per, f, time, chan, label, band):
    pow1 = psd1[(band[0] <= per) & (per <= band[1])].real
    pow2 = psd2[(band[0] <= per) & (per <= band[1])].real
    pow3 = psd3[(band[0] <= per) & (per <= band[1])].real
    f.write(label + ', 1, ' + chan + ', ' + str(time.julday).zfill(3) + ', ' + 
            str(time.hour).zfill(2) + ', ' + str(np.mean(pow1)) + ', ' + str(np.std(pow1)) + '\n')
    f.write(label + ', 2, ' + chan + ', ' + str(time.julday).zfill(3) + ', ' + 
            str(time.hour).zfill(2) + ', ' + str(np.mean(pow2)) + ', ' + str(np.std(pow2)) + '\n')
    f.write(label + ', 3, ' + chan + ', ' + str(time.julday).zfill(3) + ', ' + 
            str(time.hour).zfill(2) + ', ' + str(np.mean(pow3)) + ', ' + str(np.std(pow3)) + '\n')
    return

def calc_azi(st, inv, f):       
    st_azi = st.copy()
    st_azi.remove_response(inv)
    st_azi.decimate(4)
    st_azi.decimate(2)
    st_azi.decimate(5)
    st_azi.detrend('constant')
    st_azi.filter('bandpass', freqmin = 1/8, freqmax=1/4)
    st_azi.taper(0.05)

    # These conventions differ from the 2017 paper
    def rot(theta,sta1, loc1, sta2, loc2):
        st1 = st_azi.select(station=sta1, location=loc1)
        st2 = st_azi.select(station=sta2, location=loc2)
        theta = theta % 360.
        cosd=np.cos(np.deg2rad(theta))
        sind=np.sin(np.deg2rad(theta))
        data1 = cosd*st1[0].data - sind*st1[1].data
        data2 = sind*st1[0].data + cosd*st1[1].data
        resi = abs(sum(data1*st2[0].data)/np.sqrt(sum(data1**2)*sum(st2[0].data**2)) -1.)
        return resi

    def rot1(theta):
        return rot(theta,'ALQ2','00', 'TST4','00')
    def rot2(theta):
        return rot(theta,'TST4','00', 'TST4','10')
    def rot3(theta):
        return rot(theta,'ALQ2','00', 'TST4','10')

    theta1 = root(rot1, 0., method = 'lm' ).x[0]
    theta2 = root(rot2, 0., method = 'lm' ).x[0]
    theta3 = root(rot3, 0., method = 'lm' ).x[0]

    f.write('Orient 1 2, ' + str(st[0].stats.starttime.julday).zfill(3) + ', '
        + str(st[0].stats.starttime.hour).zfill(2) + ', ' + str(theta1) + ', '
        + str(rot1(theta1)) + '\n')
    f.write('Orient 2 3, ' + str(st[0].stats.starttime.julday).zfill(3) + ', '
        + str(st[0].stats.starttime.hour).zfill(2) + ', ' + str(theta2) + ', '
        + str(rot2(theta2)) + '\n')
    f.write('Orient 1 3, ' + str(st[0].stats.starttime.julday).zfill(3) + ', '
        + str(st[0].stats.starttime.hour).zfill(2) + ', ' + str(theta3) + ', '
        + str(rot3(theta3)) + '\n')
    return

resps = glob.glob('metadata/RESP*')
for resp in resps:
    if 'inv' not in vars():
        inv = read_inventory(resp)
    else:
        inv += read_inventory(resp)


for stime, etime in zip(stimes, etimes):
    # This loop represents a single test (we will parse data in windows for each test)
    ctime = stime
    st = Stream()
    while ctime <= etime:
        st += read('/msd/XX_TST4/2020/' + str(ctime.julday).zfill(3) + '/*BH*')
        st += read('/msd/GS_ALQ2/2020/' + str(ctime.julday).zfill(3) + '/00_BH*')
        ctime += 24*60*60
    st.merge(fill_value=0)
    st.trim(stime, etime)
    st.sort()
    if debug:
        print(st)
    f = open('Results_' + str(stime.julday).zfill(3) + '.csv','w')

    # For each data window we want to calculate the self-noise, psd, orientation
    for st_wind in st.slide(windows, windows):
        if debug:
            print(st_wind)
        for chan in ['Z', '1','2']:
            st_chan = st_wind.select(component=chan)

            p11, _ = cp(st_chan[0], st_chan[0])
            p22, _ = cp(st_chan[1], st_chan[1])
            p33, _ = cp(st_chan[2], st_chan[2])

            p21, _ = cp(st_chan[1], st_chan[0])
            p13, _ = cp(st_chan[0], st_chan[2])
            p23, fre = cp(st_chan[1], st_chan[2])
            per = 1./fre
            for idx, tr in enumerate(st_chan):
                resp = inv.get_response(tr.id, stime)
                tf, _ = resp.get_evalresp_response(1/tr.stats.sampling_rate, length, output='ACC')
                tf = tf[1:]
                
                
                if idx == 0:
                    n11 = (p11 - p21*p13/p23)/np.abs(tf)**2
                    n11 = 10*np.log10(n11)
                    psd1 = 10*np.log10(p11/np.abs(tf)**2)
                elif idx == 1:
                    n22 = (p22 - np.conjugate(p23)*p21/np.conjugate(p13))/np.abs(tf)**2
                    n22 = 10*np.log10(n22)
                    psd2 = 10*np.log10(p22/np.abs(tf)**2)
                else:
                    n33 = (p33 - p23*np.conjugate(p13)/p21)/np.abs(tf)**2
                    n33 = 10*np.log10(n33)
                    psd3 = 10*np.log10(p33/np.abs(tf)**2)
            if plot:
                plot_psd_noise(psd1, n11, psd2, n22, psd3, n33, per, st_chan)
                # Write a routine for plotting
            # Put this in a loop to make bands
            write_ratios(psd1, psd2, psd3,per, f, st_chan[0].stats.starttime, chan)
            label = 'Power .1 to 30'
            band = [0.1, 30.]
            write_noise(psd1, psd2, psd3,per, f, st_chan[0].stats.starttime, chan, label, band)
            label = 'Self-Noise .1 to 30'
            band = [0.1, 30.]
            write_noise(n11, n22, n33, per, f, st_chan[0].stats.starttime, chan, label, band)

        # Now we want to estimate the orientations for this time period
        calc_azi(st_wind, inv, f)

    f.close()
