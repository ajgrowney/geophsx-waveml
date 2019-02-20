import sys
import warnings
import glob
from obspy.clients.fdsn import Client
from obspy.core.event import Catalog
from obspy import read_events
from obspy import UTCDateTime, Stream, read
from multiprocessing import cpu_count
from eqcorrscan.utils.catalog_utils import filter_picks

from eqcorrscan.core import template_gen, match_filter, lag_calc
from eqcorrscan.utils import pre_processing, catalog_utils, plotting


def mktemplates(network_code='GEONET',
                publicIDs=['2016p008122', '2016p008353', '2016p008155',
                           '2016p008194'], plot=True):
    """Functional wrapper to make templates"""
    # We want to download some QuakeML files from the New Zealand GeoNet
    # network, GeoNet currently doesn't support FDSN event queries, so we
    # have to work around to download quakeml from their quakeml.geonet site.

    client = Client(network_code)
    # We want to download a few events from an earthquake sequence, these are
    # identified by publiID numbers, given as arguments

    catalog = Catalog()
    for publicID in publicIDs:
        if network_code == 'GEONET':
            data_stream = client._download(
                'http://quakeml.geonet.org.nz/quakeml/1.2/' + publicID)
            data_stream.seek(0, 0)
            catalog += read_events(data_stream, format="quakeml")
            data_stream.close()
        else:
            catalog += client.get_events(
                eventid=publicID, includearrivals=True)

    # Lets plot the catalog to see what we have
    if plot:
        catalog.plot(projection='local', resolution='h')

    # We don't need all the picks, lets take the information from the
    # five most used stations - note that this is done to reduce computational
    # costs.
    catalog = filter_picks(catalog, top_n_picks=5)
    # We only want the P picks in this example, but you can use others or all
    #  picks if you want.
    for event in catalog:
        for pick in event.picks:
            if pick.phase_hint == 'S':
                event.picks.remove(pick)

    # Now we can generate the templates
    templates = template_gen.template_gen(
        method='from_client', catalog=catalog, client_id=network_code,
        lowcut=2.0, highcut=9.0, samp_rate=20.0, filt_order=4, length=3.0,
        prepick=0.15, swin='all', process_len=3600, debug=0, plot=plot)

    # We now have a series of templates! Using Obspy's Stream.write() method we
    # can save these to disk for later use.  We will do that now for use in the
    # following tutorials.
    for i, template in enumerate(templates):
        template.write('tutorial_template_' + str(i) + '.ms', format='MSEED')
        # Note that this will warn you about data types.  As we don't care
        # at the moment, whatever obspy chooses is fine.
    return


def run_matchFilter_tutorial(plot=False, process_len=3600, num_cores=cpu_count()):
    """Main function to run the tutorial dataset."""
    # First we want to load our templates
    # template_names = glob.glob('tutorial_template_*.ms')
    template_names = glob.glob('./02/*.NSN___036')
    if len(template_names) == 0:
        raise IOError('Template files not found, have you run the template ' +
                      'creation tutorial?')

    templates = [read(template_name) for template_name in template_names]
    # Work out what stations we have and get the data for them
    stations = []
    for template in templates:
        for tr in template:
            stations.append((tr.stats.station, tr.stats.channel))
    # Get a unique list of stations
    stations = list(set(stations))

    # We will loop through the data chunks at a time, these chunks can be any
    # size, in general we have used 1 day as our standard, but this can be
    # as short as five minutes (for MAD thresholds) or shorter for other
    # threshold metrics. However the chunk size should be the same as your
    # template process_len.

    # You should test different parameters!!!
    start_time = UTCDateTime(2018, 2, 24, 12, 58)
    end_time = UTCDateTime(2018, 2, 24, 12, 59)
    chunks = []
    chunk_start = start_time
    while chunk_start < end_time:
        chunk_end = chunk_start + process_len
        if chunk_end > end_time:
            chunk_end = end_time
        chunks.append((chunk_start, chunk_end))
        chunk_start += process_len

    unique_detections = []

    # Set up a client to access the GeoNet database
    client = Client("GEONET")

    # Note that these chunks do not rely on each other, and could be paralleled
    # on multiple nodes of a distributed cluster, see the SLURM tutorial for
    # an example of this.
    for t1, t2 in chunks:
        # Generate the bulk information to query the GeoNet database
        bulk_info = []
        for station in stations:
            print(station[0],station[1][0],station[1][-1],t1,t2)
            bulk_info.append(('NZ', station[0], '*',
                              station[1][0] + 'H' + station[1][-1], t1, t2))

        # TODO: Find effective method to get our waveforms in bulk
        # Note this will take a little while.
        print('Downloading seismic data, this may take a while')
        st = client.get_waveforms_bulk(bulk_info)

        # TODO: Do we need to merge the stream
        # Merge the stream, it will be downloaded in chunks
        st.merge(fill_value='interpolate')

        # Pre-process the data to set frequency band and sampling rate
        # Note that this is, and MUST BE the same as the parameters used for
        # the template creation.
        print('Processing the seismic data')
        st = pre_processing.shortproc(
            st, lowcut=2.0, highcut=9.0, filt_order=4, samp_rate=20.0,
            debug=0, num_cores=num_cores, starttime=t1, endtime=t2)
        # Convert from list to stream
        st = Stream(st)

        # Now we can conduct the matched-filter detection
        detections = match_filter.match_filter(
            template_names=template_names, template_list=templates,
            st=st, threshold=8.0, threshold_type='MAD', trig_int=6.0,
            plotvar=plot, plotdir='.', cores=num_cores, debug=0,
            plot_format='png')

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
        run_matchFilter_tutorial()
    elif sys.argv[1] == "makeTemplates":
        net_code = sys.argv[2]
        idlist = list(sys.argv)[3:]
        print(idlist)
        mktemplates()
