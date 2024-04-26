% preprocess_sleep_data.m (based on Data_core_NE_EEG.m, written by Celia and Mie,
% updated 18 January 2023 by Doug Kelley for compatibility with Unix and paths on
% BlueHive.): 
% a script for preprocessing vital signs that indicate sleep state.  
% major modification by Yue Zhao in March, 2024 to add flexibility regarding
% NE and sleep scores input. See https://github.com/yzhaoinuw/preprocess_sleep_data/tree/dev
% if you have issues accessing, please email Yue at yuezhao@rochester.edu

function [] = preprocess_sleep_data(varargin)
%% 1) Define args and params
p = inputParser;

addRequired(p, 'eeg_emg_path', @ischar);

default_fp_dir = '';
default_chan_465 = '';
default_chan_405 = '';
default_chan_ttl_pulse = '';
default_interval = [];
default_EEG_stream = '';
default_EEG_chan = nan;
default_EMG_stream = '';
default_time_correction = nan;
default_sleep_score_file = '';
default_save_path = '';
default_show_figure = false;

addParameter(p, 'fp_dir', default_fp_dir, @ischar);
addParameter(p, 'EEG_stream', default_EEG_stream, @ischar);
addParameter(p, 'EEG_chan', default_EEG_chan, @isnumeric);
addParameter(p, 'EMG_stream', default_EMG_stream, @ischar);
addParameter(p, 'chan_465', default_chan_465, @ischar);
addParameter(p, 'chan_405', default_chan_405, @ischar);
addParameter(p, 'chan_ttl_pulse', default_chan_ttl_pulse, @ischar);
addParameter(p, 'interval', default_interval, @isvector);
addParameter(p, 'time_correction', default_time_correction, @isnumeric);
addParameter(p, 'sleep_score_file', default_sleep_score_file, @ischar);
addParameter(p, 'save_path', default_save_path, @ischar);
addParameter(p, 'show_figure', default_show_figure, @islogical);

% Parse the inputs.
parse(p, varargin{:})

% Access the variables.
eeg_emg_path = p.Results.eeg_emg_path;
[~,data_name,~] = fileparts(eeg_emg_path);

fp_dir = p.Results.fp_dir;
EEG_stream = p.Results.EEG_stream;
EEG_chan = p.Results.EEG_chan;
EMG_stream = p.Results.EMG_stream;
chan_465 = p.Results.chan_465;
chan_405 = p.Results.chan_405;
chan_ttl_pulse = p.Results.chan_ttl_pulse;
interval = p.Results.interval;
time_correction = p.Results.time_correction;
sleep_score_file = p.Results.sleep_score_file;
save_path = p.Results.save_path;
show_figure = p.Results.show_figure;

% determine eeg_emg_path is .exp or from fp data  
[parent_dir, filename, ext] = fileparts(eeg_emg_path);
if strcmpi(ext, '.exp')
    exp_file_path = eeg_emg_path;
else % EEG/EMG come from FP data
    exp_file_path = '';
end

if isempty(fp_dir) 
    if isempty(exp_file_path)    
        disp('No exp or TDT files found. Please check the path provided.')
        return
    end

else % if there's fp data
    if isempty(chan_465)
        disp('please provide chan_465 for signal 465.')
        return
    end
    if isempty(chan_405)
        disp('please provide chan_405 for signal 405.')
        return
    end

    % if EMG/EEG and NE come from different sources, must provide ttl pulse
    % for synchronization
    if ~strcmp(eeg_emg_path, fp_dir) 
        if isempty(chan_ttl_pulse)
            disp('Please provide TTL pulse channel to sync EEG/EMG with the fp data')
            return
        end
    end
end    

if isempty(save_path)
    save_path = [data_name '.mat'];
end

% define the following optional variables
sleep_scores = [];
ne = [];
ne_frequency = nan;

%% 2) loading and plotting EEG and EMG raw data

% if EEG/EMG come from FP data
if isempty(exp_file_path)
    fp_data = TDTbin2mat(eeg_emg_path);
    if ~isfield(fp_data.streams, EEG_stream)
        disp('Please provide EEG stream in the TDT file.')
        return
    end
    if ~isfield(fp_data.streams, EMG_stream)
        disp('Please provide EMG stream in the TDT file.')
        return
    end

    eeg_frequency = fp_data.streams.(EEG_stream).fs; %sampling frequency for EEG signal 
    eeg_data = fp_data.streams.(EEG_stream).data; %EEG signal
    eeg = eeg_data(EEG_chan,:); %add channel (1 or 2)
    emg = fp_data.streams.(EMG_stream).data; %EMG
else
    Info=loadEXP(exp_file_path,'no');
    
    TimeReldebSec=0; %start extract data from the beginning (first bin)
    %TimeRelEndSec=inf; %inf to include all data (until last bin)
    TimeRelEndSec=Info.BinFiles.Duration; %inf to include all data (including last bin)
    
    [eeg_emg_data,time]=ExtractContinuousData([],Info,[],TimeReldebSec, TimeRelEndSec,[],1);
    
    emg = eeg_emg_data(1,1:end);
    eeg = eeg_emg_data(2,1:end);
    
    %time vector using sampling frequency
    eeg_frequency = Info.Fs;
end

eeg = single(eeg);
emg = single(emg);
time_eeg = (0:length(eeg)-1)/eeg_frequency;

%% 3) Load FP (fiber_photometry) data (batch II)

if ~isempty(fp_dir)
    if ~strcmp(fp_dir, eeg_emg_path) % if NE data come from different sources
        fp_data = TDTbin2mat(fp_dir); % data is a struct  
    end
    ne_frequency = fp_data.streams.(chan_465).fs; % sampling frequency for NE, one number
    signal_465 = fp_data.streams.(chan_465).data; % hSyn-NE, array 1-D
    signal_405 = fp_data.streams.(chan_405).data; % autofluorescence, array, 1-D
    
    % if different soruces we need to sync it with fp recording
    if ~strcmp(fp_dir, eeg_emg_path)
        % removing FP trace prior to first TTL pulse
        TTL_FP = fp_data.epocs.(chan_ttl_pulse).onset; % TTL_FP is the timestamps
        TTL_gap = diff(TTL_FP) > 5 + 1; % the interval of the pulse is 5 seconds 
        if isempty(find(TTL_gap == 1, 1))
            TTL_onset = TTL_FP(1);  % when TTL pulse train is only started once
        else 
            TTL_onset = TTL_FP(find(TTL_gap==1)+1); % when TTL pulse train is started more than once
        end
        
        first_TTL = TTL_onset(1)*ne_frequency; %sampling point # to start with
        onset_FP = round(first_TTL);
        
        signal_465 = signal_465(onset_FP:end);
        signal_405 = signal_405(onset_FP:end);
    end

    % 3b) Normalize and plot (batch II)
    MeanFilterOrder = 1000; % for smoothing
    MeanFilter = ones(MeanFilterOrder,1)/MeanFilterOrder;
    
    fs_signal = 1:1:length(signal_465);
    sec_signal = fs_signal/ne_frequency;
    
    if isempty(interval)
        reg = polyfit(signal_405(1:end), signal_465(1:end), 1);
    else
        reg = polyfit(signal_405(round(interval*ne_frequency)), signal_465(round(interval*ne_frequency)), 1);
    end
    a = reg(1);
    b = reg(2);
    controlFit = a.*signal_405 + b;
    %[p,~,mu] = polyfit(signal_405(round(mouse{5}*signal_fs));                      % for scaling and centering (matlab 2020b and later)
    %controlFit = polyval(p,signal_405,[],mu);                                      % for scaling and centering
    controlFit =  filtfilt(MeanFilter,1,double(controlFit));
    normDat = (signal_465 - controlFit)./controlFit;
    delta_465 = normDat * 100;
    
    % smoothing traces
    ne = filtfilt(MeanFilter,1,double(delta_465));
    
    % downsample NE
    ds_factor_FP = 100; % also used for plotting later (section 9b)
    ne = downsample(ne, ds_factor_FP);
    ne = single(ne);
    ds_sec_signal = downsample(sec_signal, ds_factor_FP); % for plotting
    ne_frequency = ne_frequency / ds_factor_FP;
end

%% 4) read sleep scores

% NB! If there is a systematic time lag between EEG/EMG traces and scoring adjust for it by seconds here
if ~isempty(sleep_score_file)
    % Assumption: For binary vectors index 1 = time 0-1s, index 2= time 1-2 sec, and so forth
    sleep_scores = NaN(1, ceil(numel(eeg) / eeg_frequency), 'single'); 
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
% if EEG/EMG and fp data come from different sources, align them using the
% TTL pulse.
if ~isempty(fp_dir) && ~strcmp(eeg_emg_path, fp_dir)
    % TTL pulse from FP
    TTL_pulse = eeg_emg_data(3,1:end); % the actual pulse time series
    TTL_pulse_indices = find(diff(TTL_pulse>1*10^-3)==1);
    if isempty(TTL_pulse_indices) % no ttl pulse implies no FP data recorded
        onset_EEG_ind = 1;
    else
        TTL_pulse_time = TTL_pulse_indices/eeg_frequency;
        TTL_pulse_time_diff = diff(TTL_pulse_time);
        
        TTL_pulse_time_gap = TTL_pulse_time_diff > 6;
        if isempty(find(TTL_pulse_time_gap, 1))
            onset_EEG = TTL_pulse_time(1);
        else 
            onset_EEG = TTL_pulse_time(find(TTL_pulse_time_gap, 1)+1);
        end    
        onset_EEG_ind = round(onset_EEG*eeg_frequency);
    end
    
    %Cutting EEG/EMG traces leading up to first TTL 
    % Removing first seconds of EEG and EMG raw traces to align with FP trace
    emg = emg(onset_EEG_ind:end);
    eeg = eeg(onset_EEG_ind:end);
    time_eeg = (0:length(eeg)-1)/eeg_frequency;
    
    if ~isempty(sleep_scores)
        % Remove first seconds of EEG score to align with FP trace
        sleep_scores = sleep_scores(round(onset_EEG+1):end);
    end
end

 %% 6) plot (optioinal) and save extracted data to .mat file
if show_figure
    if ~isempty(fp_dir)
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
        plot(ds_sec_signal, ne)
        title('NE2m');

        figure
        a = subplot(3,1,1);
            plot(time_eeg, emg); 
            xlabel('time (s)');
            ylabel('EMG (V)');
        b = subplot(3,1,2);
            plot(time_eeg, eeg); 
            xlabel('time (s)');
            ylabel('EEG (V)');
        c = subplot(3,1,3);
            plot(ds_sec_signal, ne); 
            xlabel('time (s)');
            ylabel('NE');
        linkaxes([a, b, c],'x');

    else
        % Plot of EEG and EMG traces
        figure
        h(1) = subplot(2,1,1);
            plot(time_eeg, emg); 
            xlabel('time (s)');
            ylabel('EMG Raw (V)');
        h(2) = subplot(2,1,2);
            plot(time_eeg, eeg); 
            xlabel('time (s)');
            ylabel('EEG Raw (V)');
        linkaxes([h(1),h(2)],'x');
    end

end

num_class = 3;
%save(save_path, "trial_eeg", "trial_emg", "trial_ne", "pred_labels", "num_class", "confidence", "-v7.3")
save(save_path, "eeg", "emg", "ne", "sleep_scores", "num_class", "eeg_frequency", "ne_frequency")
end