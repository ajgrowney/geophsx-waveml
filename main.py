import glob
import matplotlib.pyplot as plt
import numpy as np
from obspy import imaging, read

st = read('./2018-01-01-0000-00M.COMB__033')

tr = st[0]
tr_fil = tr.copy()
tr_fil.filter('bandpass', freqmin=80, freqmax=150 )

t = np.arange(0, tr.stats.npts / tr.stats.sampling_rate, tr.stats.delta)

plt.subplot(2,1,1)
plt.plot(t,tr.data, 'k')

plt.subplot(2,1,2)
plt.plot(t, tr_fil.data, 'k')

plt.show()

tr.trigger("recstalta", sta=1, lta=4)
tr.plot()