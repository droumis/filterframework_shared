**------------Documentation for Flab multi-taper spectrograms------------**

**Primary Scripts:**
1. DFScripts/dfs_riptriggeredspectrogram.m

**Primary Functions:**
1. DFAnalysis/dfa_calcriptrigspectrogram.m
2. chronux_2_00/chronux/spectral_analysis/continuous/mtspecgramtrigc.m
3. DataFilter/getconstimes.m

**Primary Data Structures**
1. ripple filtered EEG (grounded)
2. ripple times (consensus)

**Summary**
1. for each tetrode, create a full epoch spectrogram
2. for each tetrode, create ripple triggered spectrograms
3. normalize each rip-trig-spect by its full-ep-spect
4. average rip-trig-spects across tetrodes for each ripple

**Common Parameters**
1. Fs = 1500 Hz
2. frequencies = 2-350 Hz
3. 300 ms window on each side of each ripple
4. 100 ms moving window, 10 ms steps
5. time-bandwidth product (TW) = 3
6. number of tapers to use (K) = 5
7. average over events (trialave) = 0

**Output**
1. S = triggered spectrum in form of time x frequency x events
2. t = times
3. f = frequencies



