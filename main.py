import sys
import warnings
import glob
import numpy as np
import matplotlib.pyplot as plt
import obspy
from pprint import pprint
from obspy.clients.fdsn import Client
from obspy.core.event import Catalog
from obspy import read_events
from obspy import UTCDateTime, Stream, read
from multiprocessing import cpu_count
from eqcorrscan.utils.catalog_utils import filter_picks

from eqcorrscan.utils.plotting import spec_trace
from eqcorrscan.core import template_gen, match_filter, lag_calc
from eqcorrscan.utils import pre_processing, catalog_utils, plotting

def get_single_stream(path):

    st = read(path, format="MSEED")
    return (st)

def get_waveforms_bulk(folder):
    all_streams_in_folder = glob.glob('./'+folder+"/*")
    bulk_return = []
    bulk_return.append(get_single_stream(all_streams_in_folder[0]))
    bulk_return.append(get_single_stream(all_streams_in_folder[1]))
    return bulk_return

def run_matchFilter(plot=False, method="av_chan_corr", threshold=0.1, min_cc=0.5, num_cores=cpu_count()):
    """Main function to run the tutorial dataset."""

    # First we want to load our templates
    # template_names = glob.glob('tutorial_template_*.ms')
    # template_names = glob.glob('./02/*.NSN___033')
    template_names = glob.glob('./01/*.NSN___030')

    print("Template Names")
    print(template_names)
    print("")
    print("")

    if len(template_names) == 0:
        raise IOError('Template files not found, have you run the template ')

    templates = [read(template_name) for template_name in template_names]
    # Work out what stations we have and get the data for them
    stations = []
    for template in templates:
        for tr in template:
            stations.append((tr.stats.station, tr.stats.channel))
    # Get a unique list of stations
    stations = list(set(stations))

    unique_detections = []

    print("Templates")
    print(templates)
    print("")
    print("")

    # DONE: Find effective method to get our waveforms in bulk
    # Note this will take a little while.
    print('Downloading seismic data locally from 2018_01')
    streams = get_waveforms_bulk("2018_01")
    # streams = [streams[0]]

    # DONE: Do we need to merge the stream
    # Merge the stream, it will be downloaded in chunks
    for st in streams:
        st.merge(fill_value='interpolate')
	#spec_trace(st,trc='white')
    # Pre-process the data to set frequency band and sampling rate
    # Note that this is, and MUST BE the same as the parameters used for
    # the template creation.
    print('Processing the seismic data')

    print("Sampling Rates")


    picked_catalog = Catalog()
    for template, template_name in zip(templates, template_names) :
        for st in streams:
            # Now we can conduct the matched-filter detection
            st = st.select(station="WK*")    # Select specific channels
            st = st.normalize()


            # ShowPlots(template)
            # ShowPlots(st, template[0])

            # print("NEW TEMPLATE ", template)

            detections = match_filter.match_filter(
                template_names=[template_name], template_list=[template], trig_int=5.0,
                st=st, threshold=threshold, threshold_type=method, plotvar=False, cores=num_cores)


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
                    # We can plot these too
                    if plot:
                        stplot = st.copy()
                        print("\nST Detection: ",stplot)
                        template = templates[template_names.index(
                            master.template_name)]
                        print("\nTemplate Detection: ",template)
                        lags = sorted([tr.stats.starttime for tr in template])
                        maxlag = lags[-1] - lags[0]
                        stplot.trim(starttime=master.detect_time - 10,
                        endtime=master.detect_time + maxlag + 10)
                        plotting.detection_multiplot(
                            stplot, template, [master.detect_time.datetime])

                picked_catalog += lag_calc.lag_calc(
                            detections=unique_detections, detect_data=st,
                            template_names=[template_name], templates=[template],
                            shift_len=1.0, min_cc=min_cc, interpolate=False, plot=False,
                            parallel=True, debug=3)
                print(len(picked_catalog))
    print('We made a total of ' + str(len(unique_detections)) + ' detections')

    return unique_detections, picked_catalog

def ShowTemplates():
    #steams and templates
    print("FINISHING TEMPLATE PLOTS")

def ShowPlots(stream, template):
    for temp in template:
        tr_filt = temp.copy()

        tr_filt.filter('lowpass', freq=1.0, corners=2, zerophase=True)

        t = np.arange(0, temp.stats.npts / temp.stats.sampling_rate, temp.stats.delta)
        plt.subplot(211)
        plt.plot(t, tr_filt.data, 'k')
        plt.ylabel('Template')
        plt.xlabel('Time [s]')

        plt.subplot(212)
        for tr in stream:
            tr_filt_st = tr.copy()

            tr_filt_st.filter('lowpass', freq=1.0, corners=2, zerophase=True)

            t_s = np.arange(0, tr.stats.npts / tr.stats.sampling_rate, tr.stats.delta)
            plt.plot(t_s, tr_filt_st.data, 'k')
            plt.ylabel('Lowpassed Stream Data')
            plt.xlabel('Time [s]')
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

        detections, picks = run_matchFilter(method=method, threshold=threshold, min_cc=0.5)
        analyzeDetections(detections)
    elif sys.argv[1] == "bulk":
        get_waveforms_bulk("2018_01")
