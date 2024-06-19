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
    'fp_dir', 'C:\Users\yzhao\matlab_projects\sleep_data_extraction\opto_NOR\408_yfp', ...
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

The first parameter is the path to your EEG/EMG data. If your EEG/EMG data was recorded with an .exp file,  then this should be the path to the .exp file here. If your EEG/EMG data was recorded using TDT system, then this should be the **directory** of the TDT files here.

The rest of the parameters are optional name-value pairs. 'interval' and 'time_correction' are always optional. The last parameter, 'show_figure' is for visualizing the preprocessed data. If you want to see it, set it to 'true'. If you don't need it, just skip it. You may or may not need to provide the other optional parameters depending on your data source. The order in which you write them does not matter. See the flowchart on the next page for the parameters you need.
![](illustration.png "Flowchart")

#### Output
The following fields are saved in the .mat file
1. eeg
2. emg
3. ne
4. sleep_scores 
5. num_class (fixed to be 3 for now)
6. eeg_frequency 
7. ne_frequency

#### Note
1. *ne* is downsampled by a factor of 100 before being saved.
2. eeg, emg, and ne are 1D arrays.
3. Make sure that your fp_dir only contains one set of TDT files. 

#### More Examples
1. EEG/EMG from .exp file, no NE data
```matlab
preprocess_sleep_data( ...
    'C:\Users\yzhao\matlab_projects\sleep_data_extraction\opto_NOR\408_yfp\408_YFP_NOR\408_YFP_NOR_2021-03-08_10-29-08-223.exp', ...
);

```

2. EEG/EMG from .exp file + NE also from .exp files.
```matlab
preprocess_sleep_data( ...
    'C:\Users\yzhao\matlab_projects\sleep_data_extraction\opto_NOR\408_yfp\408_YFP_NOR\408_YFP_NOR_2021-03-08_10-29-08-223.exp', ...
    'fp_dir', 'C:\Users\yzhao\matlab_projects\sleep_data_extraction\opto_NOR\408_yfp', ...
    'chan_465', 'BuBo', ...
    'chan_405', 'V_Bo', ...
    'chan_ttl_pulse', 'Pu2_', ...
    'interval', (100:7000), ...
    'time_correction', 2, ...
    'save_path', '408_yfp.mat', ...
    'show_figure', true ...
);

```

3. EEG/EMG from .exp file + NE from TDT files.
```matlab
preprocess_sleep_data(...
    'C:\Users\yzhao\matlab_projects\sleep_data_extraction\Day 3 sleepscore', ...
    'fp_dir', 'C:\Users\yzhao\matlab_projects\sleep_data_extraction\Day 3 sleepscore', ...
    'EEG_stream', 'EEGw', ...
    'EEG_chan', 1, ...
    'EMG_stream', 'EMG1', ...
    'chan_465', 'x465A', ...
    'chan_405', 'x405A', ...
    'interval', (1000:11000), ...
    'time_correction', 0, ...
    'show_figure', true, ...
    'save_path', 'C:\Users\yzhao\matlab_projects\sleep_data_extraction\Klaudia_709_Day3.mat' ...
);

```

4. if you already have the sleep scores (this happens if you want to contribute training data for improving the sleep scoring models).
```matlab
preprocess_sleep_data( ...
    'C:\Users\yzhao\matlab_projects\sleep_data_extraction\20221221_adra_1_238_2_242', ...
    'fp_dir', 'C:\Users\yzhao\matlab_projects\sleep_data_extraction\20221221_adra_1_238_2_242', ...
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
