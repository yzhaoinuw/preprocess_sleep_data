## preprocess_sleep_data

#### Usage example
First, make sure you "import" the preprocess_sleep_data folder and its helper, the EEGtoolbox folder, if the preprocess_sleep_data folder is not in your current directory. You could use absolute path in addpath. For example, 
```matlab
addpath('C:\Users\yzhao\matlab_projects\preprocess_sleep_data')
addpath('C:\Users\yzhao\matlab_projects\preprocess_sleep_data\EEGtoolbox')
```

After you have taken care of the "import", you can call the function like this 
```matlab
preprocess_sleep_data( ...
    'C:\Users\yzhao\matlab_projects\sleep_data_extraction\opto_NOR\408_yfp\408_YFP_NOR\408_YFP_NOR_2021-03-08_10-29-08-223.exp', ...
    'ne_dir', 'C:\Users\yzhao\matlab_projects\sleep_data_extraction\opto_NOR\408_yfp', ...
    'chan_465', 'BuBo', ...
    'chan_405', 'V_Bo', ...
    'chan_ttl_pulse', 'Pu2_', ...
    'interval', (100:7000), ...
    'time_correction', 2, ...
    'sleep_score_file', 'C:\Users\yzhao\matlab_projects\sleep_data_extraction\opto_NOR\408_yfp\408_YFP_NOR\408_NOR_sleepscore.xlsx', ...
    'save_path', '408_yfp.mat', ...
    'show_figure', true ...
);

```

The first parameter is the path to your EEG/EMG data. If your EEG/EMG data was recorded with Viewpoint,  then this should be the path to the .exp file. If your EEG/EMG data was recorded using the TDT system, then this should be the directory of the TDT files here.

The rest of the parameters are optional name-value pairs. **'interval'** and **'time_correction'** are always optional. The last parameter, **'show_figure'** is for visualizing the preprocessed data. If you want to see it, set it to **'true'**. If you don't need it, just skip it. You may or may not need to provide the other optional parameters depending on your data source. The order in which you write them does not matter. See the flowchart below for the parameters you need. Also check Section [More Examples](#More-Examples) for example usage.
![](illustration.png "Flowchart")

#### Output
The following fields are saved in the .mat file
1. eeg
2. emg
3. ne
4. sleep_scores 
5. start_time (0 if only one .bin file or for 12-hour recordings or less. See [Note](#Note) below)
6. video_start_time (TTL pulse onset in EEG signal)
7. num_class (fixed to be 3 for now)
8. eeg_frequency 
9. ne_frequency
10. video_name
11. video_path (path to .avi file if found in .exp file)


#### Note
1. *ne* is downsampled by a factor of 100 before being saved.
2. *eeg*, *emg*, and *ne* are 1D arrays.
3. Make sure that your *ne_dir* only contains one set of TDT files. 
4. Long recordings. When you have EEG/EMG data from Viewpoint and you have multiple bin files, the data in each bin file will be written to a .mat file named after the bin file. When you have EEG/EMG from TDT and the recording is longer than 12 hours, your recording will be segmented into 12-hour bins and the resulting .mat files will be named with a suffix accordingly. The **_start_time_** will mark the start time of the segmented .mat file.

#### More Examples
1. EEG/EMG from .exp file, no NE data
```matlab
preprocess_sleep_data( ...
    'C:\Users\yzhao\matlab_projects\sleep_data_extraction\opto_NOR\408_yfp\408_YFP_NOR\408_YFP_NOR_2021-03-08_10-29-08-223.exp', ...
);

```

2. EEG/EMG from .exp file + NE from TDT files.
```matlab
preprocess_sleep_data( ...
    'C:\Users\yzhao\matlab_projects\sleep_data_extraction\opto_NOR\408_yfp\408_YFP_NOR\408_YFP_NOR_2021-03-08_10-29-08-223.exp', ...
    'ne_dir', 'C:\Users\yzhao\matlab_projects\sleep_data_extraction\opto_NOR\408_yfp', ...
    'chan_465', 'BuBo', ...
    'chan_405', 'V_Bo', ...
    'chan_ttl_pulse', 'Pu2_', ...
    'save_path', '408_yfp.mat', ...
    'show_figure', true ...
);

```

3. EEG/EMG + NE from TDT files.
```matlab
preprocess_sleep_data( ...
    'C:\Users\yzhao\matlab_projects\sleep_data_extraction\830', ...
    'ne_dir', 'C:\Users\yzhao\matlab_projects\sleep_data_extraction\830', ...
    'EEG_stream', 'EEGw', ...
    'EEG_chan', 1, ...
    'EMG_stream', 'EMG1', ...
    'chan_465', 'x465C', ...
    'chan_405', 'x405C', ...
    'save_path', 'C:\Users\yzhao\matlab_projects\sleep_data_extraction\830.mat', ...
    'show_figure', true ...
);

```

4. EEG/EMG from TDT files, no NE.
```matlab
preprocess_sleep_data( ...
    'C:\Users\yzhao\matlab_projects\sleep_data_extraction\20241113_1_263_2_259_24h_test', ...
    'EEG_stream', 'EEG1', ...
    'EEG_chan', 1, ...
    'EMG_stream', 'EMG1', ...
    'save_path', 'C:\Users\yzhao\matlab_projects\sleep_data_extraction\20241113_1_263_2_259_24h_test.mat', ...
    'show_figure', true ...
);


```

5. if you already have the sleep scores (this happens if you want to contribute training data for improving the sleep scoring models).
```matlab
preprocess_sleep_data( ...
    'C:\Users\yzhao\matlab_projects\sleep_data_extraction\20221221_adra_1_238_2_242', ...
    'ne_dir', 'C:\Users\yzhao\matlab_projects\sleep_data_extraction\20221221_adra_1_238_2_242', ...
    'EEG_stream', 'EEGw', ...
    'EEG_chan', 1, ...
    'EMG_stream', 'EMG1', ...
    'chan_465', 'x465A', ...
    'chan_405', 'x405A', ...
    'interval', (1:8000), ...
    'time_correction', 0, ...
    'sleep_score_file', 'C:\Users\yzhao\matlab_projects\sleep_data_extraction\20221221_adra_1_238_2_242\211222_238_adjusted_score.xlsx', ...
    'save_path', 'C:\Users\yzhao\matlab_projects\sleep_data_extraction\20221221_adra_1_238_2_242.mat', ...
    'show_figure', true ...
);
```
