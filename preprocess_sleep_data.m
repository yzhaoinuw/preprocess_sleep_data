% preprocess_sleep_data.m (based on Data_core_NE_EEG.m, written by Celia and Mie,
% updated 18 January 2023 by Doug Kelley for compatibility with Unix and paths on
% BlueHive.): 
% a script for preprocessing vital signs that indicate sleep state.  
% major modification by Yue in March, 2024 to add flexibility regarding
% NE and sleep scores input. See https://github.com/yzhaoinuw/preprocess_sleep_data/tree/dev

function [] = preprocess_sleep_data(varargin)
p = inputParser;
%% 1) Define mouse data
% required arguments:
    % 1) EEG & EMG data folder
    % 2) 465 channel name
    % 3) 405 channel name
    % 4) TTL channel name
    % 5) Interval used for signal normalization (only batch II)
addRequired(p, 'EEG_EMG_dir', @ischar);
addRequired(p, 'chan_A', @ischar);
addRequired(p, 'chan_B', @ischar);
addRequired(p, 'chan_C', @ischar);
addRequired(p, 'interval');

default_EEG_sampling_rate = 512;
default_EMG_sampling_rate = 512;
default_NE_sampling_rate = 1017;
default_save_path = 'data.mat';
default_time_correction = 0;
default_show_figure = false;

% optional args:
    % 1) EEG sampling rate
    % 2) EMG sampling rate
    % 3) NE sampling rate
    % 4) time correction in seconds
    % 5) save path

addParameter(p, 'EEG_sampling_rate', default_EEG_sampling_rate, @isnumeric);
addParameter(p, 'EMG_sampling_rate', default_EMG_sampling_rate, @isnumeric);
addParameter(p, 'NE_sampling_rate', default_NE_sampling_rate, @isnumeric);
addParameter(p, 'time_correction', default_time_correction, @isnumeric);
addParameter(p, 'save_path', default_save_path, @ischar);
addParameter(p, 'show_figure', default_show_figure, @islogical);

% Parse the inputs.
parse(p, varargin{:})

% Access the variables.
EEG_EMG_dir = p.Results.EEG_EMG_dir;
chan_A = p.Results.chan_A;
chan_B = p.Results.chan_B;
chan_C = p.Results.chan_C;
interval = p.Results.interval;

EEG_sampling_rate = p.Results.EEG_sampling_rate;
EMG_sampling_rate = p.Results.EMG_sampling_rate;
NE_sampling_rate = p.Results.NE_sampling_rate;
time_correction = p.Results.time_correction;
save_path = p.Results.save_path;
show_figure = p.Results.show_figure;

FP_data_dir = '';
sleep_score_file = '';

% automatically look for NE data file and sleep scores
[parent_dir, ~, ~] = fileparts(EEG_EMG_dir);
parent_dir_file_list = dir(fullfile(parent_dir, '**', '*.*'));  %get list of files and folders in any subfolder
parent_dir_file_list = parent_dir_file_list(~[parent_dir_file_list.isdir]);  %remove folders from list
parent_dir_file_names = {parent_dir_file_list.name};

is_tev = endsWith(parent_dir_file_names, '.tev'); % Logical array indicating files that end with .txt
tev_file = parent_dir_file_list(is_tev);
if ~isempty(tev_file)
    FP_data_dir = tev_file.folder;
else
    disp('No NE data found.')
end

EEG_EMG_dir_file_list = dir(EEG_EMG_dir);
EEG_EMG_dir_file_list = EEG_EMG_dir_file_list(~[EEG_EMG_dir_file_list.isdir]);
EEG_EMG_dir_file_names = {EEG_EMG_dir_file_list.name};
is_exp = endsWith(EEG_EMG_dir_file_names, '.exp'); % Logical array indicating files that end with .txt
exp_file = EEG_EMG_dir_file_list(is_exp);
EEG_EMG_data_path = fullfile(exp_file.folder, exp_file.name);

is_sleepscore = contains(EEG_EMG_dir_file_names, 'score') & endsWith(EEG_EMG_dir_file_names, '.xlsx'); % Logical array indicating files that end with .txt
sleepscore = EEG_EMG_dir_file_list(is_sleepscore);
if ~isempty(sleepscore)
    sleep_score_file = fullfile(sleepscore.folder, sleepscore.name);
else
    disp('No sleep scores found.')
end

% define the following optional variables
trial_ne = [];
sleep_scores = [];
nrows_ne = Inf;

%% 2b) Load FP (fiber_photometry) data (batch II)

if ~isempty(FP_data_dir)
    data = TDTbin2mat(FP_data_dir); % data is a struct
    signal_fs = data.streams.(chan_A).fs; % sampling frequency for NE, one number
    signal_465 = data.streams.(chan_A).data; % hSyn-NE, array 1-D
    signal_405 = data.streams.(chan_B).data; % autofluorescence, array, 1-D
    
    % removing FP trace prior to first TTL pulse
    TTL_FP = data.epocs.(chan_C).onset; % TTL_FP is the timestamps
    TTL_gap = diff(TTL_FP) > 5 + 1; % the interval of the pulse is 5 seconds 
    if isempty(find(TTL_gap == 1, 1))
        TTL_onset = TTL_FP(1);  % when TTL pulse train is only started once
    else 
        TTL_onset = TTL_FP(find(TTL_gap==1)+1); % when TTL pulse train is started more than once
    end
    
    first_TTL = TTL_onset(1)*signal_fs; %sampling point # to start with
    onset_FP = round(first_TTL);
    
    signal_465 = signal_465(onset_FP:end);
    signal_405 = signal_405(onset_FP:end);

    % 3b) Normalize and plot (batch II)
    MeanFilterOrder = 1000; % for smoothing
    MeanFilter = ones(MeanFilterOrder,1)/MeanFilterOrder;
    
    fs_signal = 1:1:length(signal_465);
    sec_signal = fs_signal/signal_fs;
    
    reg = polyfit(signal_405(round(interval*signal_fs)), signal_465(round(interval*signal_fs)), 1);
    a = reg(1);
    b = reg(2);
    controlFit = a.*signal_405 + b;
    %[p,~,mu] = polyfit(signal_405(round(mouse{5}*signal_fs));                      % for scaling and centering (matlab 2020b and later)
    %controlFit = polyval(p,signal_405,[],mu);                                      % for scaling and centering
    controlFit =  filtfilt(MeanFilter,1,double(controlFit));
    normDat = (signal_465 - controlFit)./controlFit;
    delta_465 = normDat * 100;
    
    % smoothing traces
    delta465_filt = filtfilt(MeanFilter,1,double(delta_465));
    
    % downsampling traces for plotting
    ds_factor_FP = 100; % also used for plotting later (section 9b)
    ds_delta465_filt = downsample(delta465_filt, ds_factor_FP);
    ds_sec_signal = downsample(sec_signal, ds_factor_FP); % for plotting

    nrows_ne = fix(numel(delta465_filt)/NE_sampling_rate);

end

%% 3) loading and plotting EEG and EMG raw data

% Import EEG raw data to matlab
Info=loadEXP(EEG_EMG_data_path,'no');

TimeReldebSec=0; %start extract data from the beginning (first bin)
%TimeRelEndSec=inf; %inf to include all data (until last bin)
TimeRelEndSec=Info.BinFiles.Duration; %inf to include all data (including last bin)

[Data,Time]=ExtractContinuousData([],Info,[],TimeReldebSec, TimeRelEndSec,[],1);

EMG_rawtrace = Data(1,1:end);
EEG_rawtrace = Data(2,1:end);

%time vector using sampling frequency
sampling_freq = Info.Fs;
EEG_time = (0:length(EEG_rawtrace)-1)/sampling_freq;

%% 4) read sleep scores

% NB! If there is a systematic time lag between EEG/EMG traces and scoring adjust for it by seconds here
if ~isempty(sleep_score_file)
    % Assumption: For binary vectors index 1 = time 0-1s, index 2= time 1-2 sec, and so forth
    sleep_scores = zeros(1, int64(numel(EEG_rawtrace) / Info.Fs)); 
    EEG_sleepscore = readmatrix(sleep_score_file); % xlsread is not recommended by matlab, using readmatrix instead
    
    % Create binary vectors for sleep stages
    %Awake
    wake_scores = rmmissing(EEG_sleepscore(:, [2 3]));
    wake_onset = wake_scores(:, 1);
    wake_duration = wake_scores(:, 2);
    
    %Slow-wave sleep
    sws_scores = rmmissing(EEG_sleepscore(:, [6 7]));
    sws_onset = sws_scores(:, 1);
    sws_duration = sws_scores(:, 2);
    
    %REM
    rem_scores = rmmissing(EEG_sleepscore(:, [10 11]));
    REM_onset = rem_scores(:, 1);
    REM_duration = rem_scores(:, 2);

    % Most EEG scorings don't start at time 0 - which shifts the timeline of the
    % scoring relative to the EEG/EMG traces - this is corrected for below
    EEG_scoring_onset =  min([wake_onset(1), sws_onset(1)]); % the sleep scoring software shifts the sleep scores by a variable amount

    if ~isempty(REM_onset)
        EEG_scoring_onset = min([EEG_scoring_onset, REM_onset(1)]);
        REM_onset = REM_onset - EEG_scoring_onset;
        REM_onset = REM_onset + time_correction; 
        for i=1:length(REM_onset) % making time vector for EEG scoring (frequency = 1Hz)
            t = REM_onset(i)+1;
            d = REM_duration(i)-1;
            sleep_scores(t:t+d) = 2;
        end
    end

    wake_onset = wake_onset - EEG_scoring_onset;
    sws_onset = sws_onset - EEG_scoring_onset;
    
    % NB! all EEG/EMG traces are not aligned properly with sleep score (~4 s delay)
    wake_onset = wake_onset + time_correction; 
    sws_onset = sws_onset + time_correction; 
    
    for i=1:length(wake_onset) % making time vector for EEG scoring (frequency = 1Hz)
        t = wake_onset(i)+1; % +1 to put time 0 as index 1
        d = wake_duration(i)-1; % -1 compensates for adding 1
        sleep_scores(t:t+d) = 0;
    end
    
    for i=1:length(sws_onset) % making time vector for EEG scoring (frequency = 1Hz)
        t = sws_onset(i)+1; 
        d = sws_duration(i)-1;
        sleep_scores(t:t+d) = 1;
    end
    
end

%% 5) Alingment of EEG recording and FP recording

% TTL pulse from FP
TTL_pulse = Data(3,1:end); % the actual pulse time series
TTL_pulse_indices = find(diff(TTL_pulse>1*10^-3)==1);
if isempty(TTL_pulse_indices) % no ttl pulse implies no FP data recorded
    onset_EEG_ind = 1;
else
    TTL_pulse_time = TTL_pulse_indices/sampling_freq;
    TTL_pulse_time_diff = diff(TTL_pulse_time);
    
    TTL_pulse_time_gap = TTL_pulse_time_diff > 6;
    if isempty(find(TTL_pulse_time_gap, 1))
        onset_EEG = TTL_pulse_time(1);
    else 
        onset_EEG = TTL_pulse_time(find(TTL_pulse_time_gap, 1)+1);
    end    
    onset_EEG_ind = round(onset_EEG*sampling_freq);
end

%Cutting EEG/EMG traces leading up to first TTL 
% Removing first seconds of EEG and EMG raw traces to align with FP trace
EMG_rawtrace_cut = EMG_rawtrace(onset_EEG_ind:end);
EEG_rawtrace_cut = EEG_rawtrace(onset_EEG_ind:end);
EEG_time_cut = (0:length(EEG_rawtrace_cut)-1)/sampling_freq;

if ~isempty(sleep_scores)
    % Remove first seconds of EEG score to align with FP trace
    sleep_scores = sleep_scores(round(onset_EEG+1):end);
end

 %% 6) reshaping EEG, EMG, and NE

% divide by sampling rate to get the duration in seconds
nrows_eeg = fix(numel(EEG_rawtrace_cut)/EEG_sampling_rate); 
nrows_emg = fix(numel(EMG_rawtrace_cut)/EMG_sampling_rate);
%nrows_ne = fix(numel(delta465_filt)/NE_sampling_rate);
nrows = min([nrows_eeg nrows_emg nrows_ne]); % take the smallest duration

trial_eeg = transpose(reshape(EEG_rawtrace_cut(1:nrows*EEG_sampling_rate), [EEG_sampling_rate, nrows]));
trial_emg = transpose(reshape(EMG_rawtrace_cut(1:nrows*EMG_sampling_rate), [EMG_sampling_rate, nrows]));

if ~isempty(FP_data_dir)
    trial_ne = transpose(reshape(delta465_filt(1:nrows*NE_sampling_rate), [NE_sampling_rate, nrows]));
end

if show_figure
    if ~isempty(FP_data_dir)
        figure
        a = subplot(4,1,1);
        plot(sec_signal(1000:end), signal_405(1000:end));
        title('raw control');
        b = subplot(4,1,2);
        plot(sec_signal(1000:end), signal_465(1000:end));
        title('raw signal');
        c = subplot(4,1,3);
        plot(sec_signal(1000:end), signal_465(1000:end));
        hold on
        plot(sec_signal(1000:end), controlFit(1000:end));
        title('fitted control');
        d = subplot(4,1,4);
        plot(sec_signal(1000:end), delta_465(1000:end));
        title('normalized signal');
        linkaxes([a,b,c,d],'x');
        % Plot of the three traces above each other (the index 1000:end removes the
        % first second of the recoding for nicer plotting)
        figure
        plot(ds_sec_signal, ds_delta465_filt)
        title('NE2m');

        figure
        a = subplot(3,1,1);
            plot(EEG_time_cut, EMG_rawtrace_cut); 
            xlabel('time (s)');
            ylabel('EMG (V)');
        b = subplot(3,1,2);
            plot(EEG_time_cut, EEG_rawtrace_cut); 
            xlabel('time (s)');
            ylabel('EEG (V)');
        c = subplot(3,1,3);
            plot(ds_sec_signal, ds_delta465_filt); 
            xlabel('time (s)');
            ylabel('NE');
        linkaxes([a, b, c],'x');

    end

    % Plot of EEG and EMG traces
    figure
    h(1) = subplot(2,1,1);
        plot(EEG_time, EMG_rawtrace); 
        xlabel('time (s)');
        ylabel('EMG Raw (V)');
    h(2) = subplot(2,1,2);
        plot(EEG_time, EEG_rawtrace); 
        xlabel('time (s)');
        ylabel('EEG Raw (V)');
    linkaxes([h(1),h(2)],'x');

end

pred_labels = sleep_scores;
num_class = 3;
confidence = [];
%save(save_path, "trial_eeg", "trial_emg", "trial_ne", "sleep_scores", "-v7.3")
save(save_path, "trial_eeg", "trial_emg", "trial_ne", "pred_labels", "num_class", "confidence", "-v7.3")
end