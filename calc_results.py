#!/usr/bin/env python
import glob
import numpy as np
from obspy.core import UTCDateTime, read, Stream
from obspy.core.inventory import read_inventory
from matplotlib.mlab import csd
from scipy.optimize import root

debug, plot = True, False
stimes, etimes = [], []
length, overlap, windows = 2**12, 2**8, 3*60*60

# West vault firs times  need to adjust
stimes.append(UTCDateTime('2020-282T00:00:00'))
etimes.append(UTCDateTime('2020-283T00:00:00'))

def cp(tr1,tr2):
    cpval,fre = csd(tr1.data, tr2.data, NFFT=length,
                    Fs=tr1.stats.sampling_rate,
                    noverlap=overlap, scale_by_freq=True)
    return cpval[1:], fre[1:]

def plot_psd_noise():
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

def calc_azi(st, f):       
    st_azi = st.copy()
    st_azi.decimate(4)
    st_azi.decimate(2)
    st_azi.decimate(5)
    st_azi.detrend('constant')
    st_azi.filter('bandpass', freqmin = 1/8, freqmax=1/4)
    st_azi.taper(0.05)

    def rot(theta,sta1, loc1, sta2, loc2):
        st1 = st_azi.select(station=sta1, location=loc1)
        print(st1)
        st2 = st_azi.select(station=sta2, location=loc2)
        print(st2)
        theta = theta % 360.
        cosd=np.cos(np.deg2rad(theta))
        sind=np.sin(np.deg2rad(theta))
        data1 = cosd*st1[0].data - sind*st2[0].data
        data2 = sind*st1[0].data + cosd*st2[0].data
        resi = abs(sum(data1*st1[1].data)/np.sqrt(sum(data1**2)*sum(st1[1].data**2)) -1.)
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
        for chan in ['1','2','Z']:
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
                tf, _ = resp.get_evalresp_response(tr.stats.sampling_rate, length)
                tf = tf[1:]
                if idx == 0:
                    n11 = ((2*np.pi*fre)**2)*(p11 - p21*p13/p23)/np.abs(tf)**2
                    n11 = 10*np.log10(n11)
                    psd1 = 10*np.log10(p11/np.abs(tf)**2)
                elif idx == 1:
                    n22 = ((2*np.pi*fre)**2)*(p22 - np.conjugate(p23)*p21/np.conjugate(p13))/np.abs(tf)**2
                    n22 = 10*np.log10(n22)
                    psd2 = 10*np.log10(p22/np.abs(tf)**2)
                else:
                    n33 = ((2*np.pi*fre)**2)*(p33 - p23*np.conjugate(p13)/p21)/np.abs(tf)**2
                    n33 = 10*np.log10(n33)
                    psd3 = 10*np.log10(p33/np.abs(tf)**2)
            if plot:
                plot_psd_noise()
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
        calc_azi(st_wind, f)

    f.close()


