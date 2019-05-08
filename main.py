import sys
import warnings
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
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

from eqcorrscan.utils.plotting import spec_trace, detection_multiplot
from eqcorrscan.core import template_gen, match_filter, lag_calc
from eqcorrscan.utils import pre_processing, catalog_utils, plotting


class NullWriter(object):
    def write(self, arg):
        pass


# Error with conflicting installations of numpy
# Solution is to specify
# Solution link > https://stackoverflow.com/questions/20554074/sklearn-omp-error-15-when-fitting-models
import os
os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"

def get_single_stream(path):
    st = read(path, format="MSEED")
    return (st)

def get_waveforms_bulk(folder):
    all_streams_in_folder = glob.glob('./'+folder+"/2016-07-04*.IRIS__024")
    bulk_return = [get_single_stream(st) for st in all_streams_in_folder]
    # bulk_return.append(get_single_stream(all_streams_in_folder[0]))
    # bulk_return.append(get_single_stream(all_streams_in_folder[1]))
    return bulk_return

def filterStream(stream):

    filtered_stream = stream.filter('bandpass',freqmin=3.0, freqmax=20.0)
    return filtered_stream


def run_matchFilter(plot=False, method="av_chan_corr", threshold=0.1, min_cc=0.5, num_cores=cpu_count()):
    """Main function to run the tutorial dataset."""
    nullwrite = NullWriter()
    # First we want to load our templates
    template_names = glob.glob('./07/2016-07-04*.NSN___024')

    print("Template Names")
    print(template_names)
    print("")

    if len(template_names) == 0:
        raise IOError('Template files not found, have you run the template ')

    templates = [read(template_name) for template_name in template_names]


    # print("Templates")
    # print(templates)
    # print("")
    # DONE: Get Stream Data to compare against templates
    print('Downloading seismic data locally from 2016_07')
    streams = get_waveforms_bulk("2016_07")
    print(len(streams))
    # print("Streams")
    # print(streams[0])
    # print("...")

    # Initialize all return values
    unique_detections, picked_catalog, detectedWaveforms = [], Catalog(), []

    for template, template_name in zip(templates, template_names) :
        for st in streams:
            # Now we can conduct the matched-filter detection
            st = st.select(station="WK*")    # Select specific stations
            st = filterStream(st)

            oldout = sys.stdout
            sys.stdout = nullwrite
            detections = match_filter.match_filter(
                template_names=[template_name], template_list=[template], trig_int=5.0,
                st=st, threshold=threshold, threshold_type=method, plotvar=False, cores=num_cores)
            sys.stdout = oldout


            if len(detections) > 0:
                waveforms = match_filter.extract_from_stream(st,detections)
                detectedWaveforms += waveforms

            current_picks = Catalog()
            for master in detections:
                keep = True
                for slave in detections:
                    if not master == slave and abs(master.detect_time -
                                                   slave.detect_time) <= 1.0:
                        # If the events are within 1s of each other then test which
                        # was the 'best' match, strongest detection
                        if not master.detect_val > slave.detect_val:
                            keep = False
                            break
                if keep:
                    unique_detections.append(master)
                    print('Detection at :' + str(master.detect_time) +
                          ' for template ' + master.template_name +
                          ' with a cross-correlation sum of: ' +
                          str(master.detect_val))
                oldout = sys.stdout
                sys.stdout = nullwrite
                current_picks += lag_calc.lag_calc(
                            detections=unique_detections, detect_data=st,
                            template_names=[template_name], templates=[template],
                            shift_len=1.0, min_cc=min_cc, interpolate=False, plot=False,
                            parallel=True, debug=3)
                sys.stdout = oldout

            copy_stream, copy_template = st.copy(), template.copy()

            figures = []
            for event in current_picks:
                times = [min([pk.time -0.05 for pk in event.picks])]
                plotString = str(event.resource_id)+".png"
                mPlotFigure = detection_multiplot(stream=copy_stream, template=copy_template,times=times, size=(10.5, 7.5), savefile=plotString, return_figure=True)
                figures.append(mPlotFigure)
                print("Figure Plotted " + str(figures.count))
                picked_catalog += current_picks

            if len(current_picks) > 0:
                pdfName = str(current_picks[0].resource_id)
                pdfName.replace('.','').replace('/','')
                pdf = matplotlib.backends.backend_pdf.PdfPages(pdfName+"output.pdf")
                for fig in figures:
                    pdf.savefig( fig )
                pdf.close()

    print('We made a total of ' + str(len(unique_detections)) + ' detections')
    print('We made '+str(len(picked_catalog))+ ' picks')
    print('We have '+str(len(detectedWaveforms)) + ' detected waveforms')
    return unique_detections, picked_catalog, detectedWaveforms

def analyzeDetections(detections):
    detectValSorted = sorted(detections, key=lambda x: x.detect_val, reverse=True)
    pairs = [('Template: ' + x.template_name,'Detection: '+str(x.detect_time),'Value: '+str(x.detect_val)) for x in detections]
    [ print(pair) for pair in pairs]
if __name__ == '__main__':

    if sys.argv[1] == "matchFilter":
        method = sys.argv[2] if len(sys.argv) > 2 else 'absolute'
        threshold = float(sys.argv[3]) if len(sys.argv) > 3 else 3.0

        detections, picks,waveforms = run_matchFilter(method=method, threshold=threshold, min_cc=0.5)

        detectionData = []
        for i in range(len(detections)):
            detectionData.append([detections[i].detect_time, detections[i].template_name, method, threshold, str(detections[i].detect_val), str(len(detections)), 'Y'])
        currentPickDF = pd.DataFrame(detectionData, columns = ['time', 'template', 'method', 'thresh', 'xcorr', 'num_dets', 'class'])
        currentPickDF.to_csv("detections.csv", sep='\t', encoding='utf-8')

        analyzeDetections(detections)
    elif sys.argv[1] == "bulk":
        get_waveforms_bulk("2016_07")
