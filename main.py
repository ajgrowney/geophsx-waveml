import sys
import warnings
import glob
from pprint import pprint
from obspy.clients.fdsn import Client
from obspy.core.event import Catalog
from obspy import read_events
from obspy import UTCDateTime, Stream, read
from multiprocessing import cpu_count
from eqcorrscan.utils.catalog_utils import filter_picks

from eqcorrscan.core import template_gen, match_filter, lag_calc
from eqcorrscan.utils import pre_processing, catalog_utils, plotting

def get_single_stream(path):

    st = read(path, format="MSEED")
    # temp1 = st = Stream(st)read("./02/2018-02-14-1329-50M.NSN___036") # obspy Stream
    return (st)

def get_waveforms_bulk(folder):
    all_streams_in_folder = glob.glob('./'+folder+"/*")
    bulk_return = []
    bulk_return.append(get_single_stream(all_streams_in_folder[0]))
    return bulk_return

def run_matchFilter(plot=False, process_len=100, num_cores=cpu_count()):
    """Main function to run the tutorial dataset."""
    # First we want to load our templates
    # template_names = glob.glob('tutorial_template_*.ms')
    template_names = glob.glob('./02/*.NSN___036')
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

    # Note that these chunks do not rely on each other, and could be paralleled
    # on multiple nodes of a distributed cluster, see the SLURM tutorial for
    # an example of this.
    # for t1, t2 in chunks:
    #     # Generate the bulk information to query the GeoNet database
    #     t1Folder = str(t1.year) + "_" + (str(t1.month).zfill(2)) # '2018_01'
    #     t2Folder = str(t2.year) + "_" + (str(t2.month).zfill(2)) # '2018_01'
    #
    #     bulk_info = []
    #     for station in stations:
    #         print("Station",type(station))
    #         print("Time",type(t1),t2)
    #         print("Month t1", t1.minute, t2.minute)
    #         print(t1Folder)
    #         print(station)
    #         print(station[0])
    #         print(station[1][-1])
    #         bulk_info.append(())
    #         return

    # DONE: Find effective method to get our waveforms in bulk
    # Note this will take a little while.
    print('Downloading seismic data locally from 2018_01')
    streams = get_waveforms_bulk("2018_01")
    #streams = [streams[0]]

    # DONE: Do we need to merge the stream
    # Merge the stream, it will be downloaded in chunks
    for st in streams:
        st.merge(fill_value='interpolate')

    # Pre-process the data to set frequency band and sampling rate
    # Note that this is, and MUST BE the same as the parameters used for
    # the template creation.
    print('Processing the seismic data')

    # for st in streams:
    #     print(vars(st))
    #     return
    #     st = pre_processing.shortproc(
    #         st, lowcut=2.0, highcut=9.0, filt_order=4, samp_rate=20.0,
    #         debug=0, num_cores=num_cores, starttime=t1, endtime=t2)
    #     # Convert from list to stream
    #     st = Stream(st)

    for st in streams:
        # Now we can conduct the matched-filter detection
        print(type(st))
        st = st.select(channel="EH*")
        print(st._get_common_channels_info())
        print((st))
        print(type(templates[0]))
        [tr.decimate(factor=4, strict_length=False) for tr in st]
        templates = [st]
        template_names = ['Template 1']

        detections = match_filter.match_filter(
            template_names=template_names, template_list=templates, trig_int=1.0,
            st=st, threshold=0.5, threshold_type='MAD', plotvar=False, cores=num_cores)

        # Now lets try and work out how many unique events we have just to
        # compare with the GeoNet catalog of 20 events on this day in this
        # sequence
        for master in detections:
            keep = True
            for slave in detections:
                if not master == slave and abs(master.detect_time -
                                               slave.detect_time) <= 1.0:
                    # If the events are within 1s of each other then test which
                    # was the 'best' match, strongest detection
                    if not master.detect_val > slave.detect_val:
                        keep = False
                        print('Removed detection at %s with cccsum %s'
                              % (master.detect_time, master.detect_val))
                        print('Keeping detection at %s with cccsum %s'
                              % (slave.detect_time, slave.detect_val))
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
                    template = templates[template_names.index(
                        master.template_name)]
                    lags = sorted([tr.stats.starttime for tr in template])
                    maxlag = lags[-1] - lags[0]
                    stplot.trim(starttime=master.detect_time - 10,
                                endtime=master.detect_time + maxlag + 10)
                    plotting.detection_multiplot(
                        stplot, template, [master.detect_time.datetime])
    print('We made a total of ' + str(len(unique_detections)) + ' detections')
    return unique_detections

def run_pickAndLag_tutorial(min_magnitude=2, shift_len=0.2, num_cores=4, min_cc=0.5):
    """Functional, tested example script for running the lag-calc tutorial."""
    if num_cores > cpu_count():
        num_cores = cpu_count()
    client = Client('NCEDC')
    t1 = UTCDateTime(2004, 9, 28)
    t2 = t1 + 86400
    print('Downloading catalog')
    catalog = client.get_events(
        starttime=t1, endtime=t2, minmagnitude=min_magnitude,
        minlatitude=35.7, maxlatitude=36.1, minlongitude=-120.6,
        maxlongitude=-120.2, includearrivals=True)
    # We don't need all the picks, lets take the information from the
    # five most used stations - note that this is done to reduce computational
    # costs.
    catalog = catalog_utils.filter_picks(
        catalog, channels=['EHZ'], top_n_picks=5)
    # There is a duplicate pick in event 3 in the catalog - this has the effect
    # of reducing our detections - check it yourself.
    for pick in catalog[3].picks:
        if pick.waveform_id.station_code == 'PHOB' and \
                        pick.onset == 'emergent':
            catalog[3].picks.remove(pick)
    print('Generating templates')
    templates = template_gen.from_client(
        catalog=catalog, client_id='NCEDC', lowcut=2.0, highcut=9.0,
        samp_rate=50.0, filt_order=4, length=3.0, prepick=0.15,
        swin='all', process_len=3600)
    # In this section we generate a series of chunks of data.
    start_time = UTCDateTime(2004, 9, 28, 17)
    end_time = UTCDateTime(2004, 9, 28, 20)
    process_len = 3600
    chunks = []
    chunk_start = start_time
    while chunk_start < end_time:
        chunk_end = chunk_start + process_len
        if chunk_end > end_time:
            chunk_end = end_time
        chunks.append((chunk_start, chunk_end))
        chunk_start += process_len

    all_detections = []
    picked_catalog = Catalog()
    template_names = [str(template[0].stats.starttime)
                      for template in templates]
    for t1, t2 in chunks:
        print('Downloading and processing for start-time: %s' % t1)
        # Download and process the data
        bulk_info = [(tr.stats.network, tr.stats.station, '*',
                      tr.stats.channel, t1, t2) for tr in templates[0]]
        # Just downloading a chunk of data
        st = client.get_waveforms_bulk(bulk_info)
        st.merge(fill_value='interpolate')
        st = pre_processing.shortproc(
            st, lowcut=2.0, highcut=9.0, filt_order=4, samp_rate=50.0,
            debug=0, num_cores=num_cores)
        detections = match_filter.match_filter(
            template_names=template_names, template_list=templates, st=st,
            threshold=8.0, threshold_type='MAD', trig_int=6.0, plotvar=False,
            plotdir='.', cores=num_cores)
        # Extract unique detections from set.
        unique_detections = []
        for master in detections:
            keep = True
            for slave in detections:
                if not master == slave and\
                   abs(master.detect_time - slave.detect_time) <= 1.0:
                    # If the events are within 1s of each other then test which
                    # was the 'best' match, strongest detection
                    if not master.detect_val > slave.detect_val:
                        keep = False
                        break
            if keep:
                unique_detections.append(master)
        all_detections += unique_detections

        picked_catalog += lag_calc.lag_calc(
            detections=unique_detections, detect_data=st,
            template_names=template_names, templates=templates,
            shift_len=shift_len, min_cc=min_cc, interpolate=False, plot=False,
            parallel=True, debug=3)
    # Return all of this so that we can use this function for testing.
    return all_detections, picked_catalog, templates, template_names


if __name__ == '__main__':
    if sys.argv[1] == "pickAndLag" :
        run_pickAndLag_tutorial(min_magnitude=4, num_cores=cpu_count())
    elif sys.argv[1] == "matchFilter":
        run_matchFilter()
    elif sys.argv[1] == "bulk":
        get_waveforms_bulk("2018_01")
