import sys
import warnings
import glob
import numpy as np
import matplotlib.pyplot as plt
import re
import obspy
from pprint import pprint
from obspy.clients.fdsn import Client
from obspy.core.event import Catalog
from obspy.signal.filter import bandpass
from obspy import read_events
from obspy import UTCDateTime, Stream, read
from multiprocessing import cpu_count
from eqcorrscan.utils.catalog_utils import filter_picks

from eqcorrscan.utils.plotting import spec_trace
from eqcorrscan.core import template_gen, match_filter, lag_calc
from eqcorrscan.utils import pre_processing, catalog_utils, plotting

# Error with conflicting installations of numpy
# Solution is to specify
# Solution link > https://stackoverflow.com/questions/20554074/sklearn-omp-error-15-when-fitting-models
import os
os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"

def get_single_stream(path):

    st = read(path, format="MSEED")
    return (st)

def get_waveforms_bulk(folder):
    all_streams_in_folder = glob.glob('./'+folder+"/*")
    bulk_return = []
    bulk_return.append(get_single_stream(all_streams_in_folder[0]))
    bulk_return.append(get_single_stream(all_streams_in_folder[1]))
    return bulk_return

def filterStream(stream):

    filtered_stream = stream.filter('bandpass',freqmin=3.0, freqmax=20.0)
    return filtered_stream


def run_matchFilter(plot=False, method="av_chan_corr", threshold=0.1, min_cc=0.5, num_cores=cpu_count()):
    """Main function to run the tutorial dataset."""

    # First we want to load our templates
    template_names = glob.glob('./07/*.NSN___030')

    print("Template Names")
    print(template_names)
    print("")

    if len(template_names) == 0:
        raise IOError('Template files not found, have you run the template ')

    templates = [read(template_name) for template_name in template_names]


    print("Templates")
    print(templates)
    print("")

    # DONE: Get Stream Data to compare against templates
    print('Downloading seismic data locally from 2016_07')
    streams = get_waveforms_bulk("2016_07")
    print("Streams")
    print(streams[0])
    print("...")

    # Initialize all return values
    unique_detections, picked_catalog, detectedWaveforms = [], Catalog(), []

    for template, template_name in zip(templates, template_names) :
        for st in streams:
            # Now we can conduct the matched-filter detection
            st = st.select(station="WK*")    # Select specific stations
            st = filterStream(st)

            detections = match_filter.match_filter(
                template_names=[template_name], template_list=[template], trig_int=5.0,
                st=st, threshold=threshold, threshold_type=method, plotvar=False, cores=num_cores)

            if len(detections) > 0:
                waveforms = match_filter.extract_from_stream(st,detections)
                detectedWaveforms += waveforms
                times = [min([pk.time -0.05 for pk in event.picks])]
                detection_multiplot(stream=waveforms, template=template, times=times)
                # for w,templ in zip(waveforms,template):
                #     templ.plot()
                #     w.plot()
                    


            for master in detections:
                keep = True
                for slave in detections:
                    if not master == slave and abs(master.detect_time -
                                                   slave.detect_time) <= 1.0:
                        # If the events are within 1s of each other then test which
                        # was the 'best' match, strongest detection
                        if not master.detect_val > slave.detect_val:
                            keep = False
                            #print('Removed detection at %s with cccsum %s'
                            #      % (master.detect_time, master.detect_val))
                            #print('Keeping detection at %s with cccsum %s'
                            #      % (slave.detect_time, slave.detect_val))
                            break
                if keep:
                    unique_detections.append(master)
                    print('Detection at :' + str(master.detect_time) +
                          ' for template ' + master.template_name +
                          ' with a cross-correlation sum of: ' +
                          str(master.detect_val))
                    # # We can plot these too
                    # print("------------------------------")
                    # print(master)
                    # print(detections)
                    # tm = master.detect_time
                    # print(tm)
                    # print(tm.hour)
                    # print(tm.minute)
                    # print(tm.second)
                    # print(tm.microsecond)
                    # timeInSeconds = tm.hour*60*60+tm.minute*60+tm.second+tm.microsecond*0.000001
                    # print(timeInSeconds)
                    # print("------------------------------")
                    # count = 0
                    # for trace, trace_temp in zip(st, template):
                    #     ShowPlots(trace, trace_temp, template_name, timeInSeconds, str(len(unique_detections)), count)
                    #     count += 1
                    # if plot:
                    #     stplot = st.copy()
                    #     print("\nST Detection: ",stplot)
                    #     template = templates[template_names.index(
                    #         master.template_name)]
                    #     print("\nTemplate Detection: ",template)
                    #     lags = sorted([tr.stats.starttime for tr in template])
                    #     maxlag = lags[-1] - lags[0]
                    #     stplot.trim(starttime=master.detect_time - 10,
                    #     endtime=master.detect_time + maxlag + 10)
                    #     plotting.detection_multiplot(
                    #         stplot, template, [master.detect_time.datetime])

                picked_catalog += lag_calc.lag_calc(
                            detections=unique_detections, detect_data=st,
                            template_names=[template_name], templates=[template],
                            shift_len=1.0, min_cc=min_cc, interpolate=False, plot=False,
                            parallel=True, debug=3)
    print('We made a total of ' + str(len(unique_detections)) + ' detections')
    print('We made '+str(len(picked_catalog))+ ' picks')
    return unique_detections, picked_catalog, detectedWaveforms

def ShowPlots(stream, template, tempName, tm, detectCount, count):
    tr_filt = template.copy()

    tr_filt.filter('lowpass', freq=1.0, corners=2, zerophase=True)

    t = np.arange(0, template.stats.npts / template.stats.sampling_rate, template.stats.delta)
    ax1 = plt.subplot(211)
    # ax1.set_xlim(left=tm)

    plt.plot(t, tr_filt.data, 'k')
    plt.ylabel('Template - ' + tempName)
    plt.xlabel('Time [s]')

    ax2 = plt.subplot(212)
    tr_filt_st = stream.copy()

    tr_filt_st.filter('lowpass', freq=5.0, corners=2, zerophase=True)
    ax2.set_xlim(left=tm, right=(tm+160))

    t_s = np.arange(0, tr_filt_st.stats.npts / tr_filt_st.stats.sampling_rate, tr_filt_st.stats.delta)
    plt.plot(t_s, tr_filt_st.data, 'k')
    plt.ylabel('Lowpassed Stream Data')
    plt.xlabel('Time [s]')
    plt.savefig('./plots/'+detectCount+'-'+str(count)+'.jpg')
    plt.show()
    print("FINISHING STREAM PLOTS")

def analyzeDetections(detections):
    detectValSorted = sorted(detections, key=lambda x: x.detect_val, reverse=True)
    pairs = [('Template: ' + x.template_name,'Detection: '+str(x.detect_time),'Value: '+str(x.detect_val)) for x in detections]
    [ print(pair) for pair in pairs]
if __name__ == '__main__':

    if sys.argv[1] == "matchFilter":
        method = sys.argv[2] if len(sys.argv) > 2 else 'absolute'
        threshold = float(sys.argv[3]) if len(sys.argv) > 3 else 3.0

        detections, picks,waveforms = run_matchFilter(method=method, threshold=2.0, min_cc=0.5)

        analyzeDetections(detections)
    elif sys.argv[1] == "bulk":
        get_waveforms_bulk("2016_07")
