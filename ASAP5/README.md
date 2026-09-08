# ASAP5 voltage sensor imaging in *ex vivo* mouse brain slices

**Installation:**

To install MATLAB, you require a paid MATLAB license or subscription for both PC (Windows 10 or 11) and MacOS machines. The MATLAB code in this folder has been tested on PC machines (Windows 10 or 11) using MATLAB R2024a. You can install Python [here](https://www.python.org/downloads/). The Python code in this folder as been tested on MacOS machine.

**Files:**

The following files are included in this folder:

Voltage_sensor_analysis.m: Script to process (correct for bleaching, flip and extract dF/F0) of fluorescence changes of the ASAP5 sensor in response to electrical stimulations.

Voltage_sensor_analysis_for_2p_different_IPI.m: Script to process (correct for bleaching, flip and extract dF/F0) of fluorescence changes of the ASAP5 sensor in response to 2 pulses of electrical stimulations at different interpuls intervals.

Voltage_sensor_spontaneous_nAChR_event_detection.ipynb: Script to detect and extract information (amplitude, duration, etc...) about spontaneous nAChR events 

Voltage_Sensor_analysis_nAChR_mediated_depolarisation_metrics.ipynb: Script to extract metrics of nAChR-mediated depolarisation (extracted by subtraction of responses to 1 electrical pulse in and out nAChR blocker DHBE)

**Instructions:**

Before using the the Voltage_sensor_analysis.m and Voltage_sensor_analysis_for_2p_different_IPI.m codes you need to create, in MATLAB, two tables as following:

1) an input table called A with: Column 1 = with the frame number Column 2 until the end of experiment = raw fluorescence transients extracted from imageJ (1 column per recording).

AND

2) an input table called "stims" (1 line only) with: Column 1 = 0 Column 2 until the end of experiment= the stimulation done for the corresponding fluoresence transients (1 column per recording).


For the Voltage_sensor_spontaneous_nAChR_event_detection.ipynb and Voltage_Sensor_analysis_nAChR_mediated_depolarisation_metrics.ipynb, data will be extracted directly from the raw excel files organised with culumn headers as recording number & electrical stimulation applied for that recording, with below the raw fluorescence extracted from ImageJ.
