jGCaMP8m sensor imaging in ex vivo mouse brain slices

**Installation:**

To install MATLAB on a PC you require a paid subscription from MATLAB.


**Files:**

The following file is included in this folder:

jGCaMP8m_analysis.m: Script to process (correct for bleaching and extract dF/F0) of fluorescence changes of the jGCaMP8m sensor in response to electrical stimulations.


**Instructions:**

Before using the code you need to create two tables as following:

an input table called A with:
Column 1 = with the frame number 
Column 2 until the end of experiment = raw fluorescence transients extracted from imageJ (1 column per recording)

AND

an input table called "stims" (1 line only) with:
Column 1 = 0 
Column 2 until the end of experiment= the stimulation done for the corresponding fluoresence transients (1 column per recording)
