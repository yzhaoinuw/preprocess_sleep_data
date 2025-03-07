% preprocess_sleep_data.m (based on Data_core_NE_EEG.m, written by Celia and Mie,
% updated 18 January 2023 by Doug Kelley for compatibility with Unix and paths on
% BlueHive.): 
% a script for preprocessing vital signs that indicate sleep state.  
% major modification by Yue Zhao since March, 2024 to add flexibility regarding
% NE and sleep scores input. See https://github.com/yzhaoinuw/preprocess_sleep_data/tree/dev
% if you have issues accessing, please email Yue at yuezhao@rochester.edu
% $Version: 0.2.6

function [] = preprocess_sleep_data(varargin)
%% 1) Define args and params
p = inputParser;

addRequired(p, 'eeg_emg_path', @ischar);

default_ne_dir = '';
default_chan_465 = '';
default_chan_405 = '';
default_chan_ttl_pulse = '';
default_interval = [];
default_EEG_stream = '';
default_EEG_chan = 1;
default_EMG_stream = '';
default_ds_factor_FP = 100;
default_time_correction = 0;
default_sleep_score_file = '';
default_save_path = '';
default_show_figure = false;

addParameter(p, 'ne_dir', default_ne_dir, @ischar);
addParameter(p, 'EEG_stream', default_EEG_stream, @ischar);
addParameter(p, 'EEG_chan', default_EEG_chan, @isnumeric);
addParameter(p, 'EMG_stream', default_EMG_stream, @ischar);
addParameter(p, 'chan_465', default_chan_465, @ischar);
addParameter(p, 'chan_405', default_chan_405, @ischar);
addParameter(p, 'chan_ttl_pulse', default_chan_ttl_pulse, @ischar);
addParameter(p, 'interval', default_interval, @isvector);
addParameter(p, 'ds_factor_FP', default_ds_factor_FP, @isnumeric);
addParameter(p, 'time_correction', default_time_correction, @isnumeric);
addParameter(p, 'sleep_score_file', default_sleep_score_file, @ischar);
addParameter(p, 'save_path', default_save_path, @ischar);
addParameter(p, 'show_figure', default_show_figure, @islogical);

% Parse the inputs.
parse(p, varargin{:})

% Access the variables.
eeg_emg_path = p.Results.eeg_emg_path;
[~,data_name,~] = fileparts(eeg_emg_path);

ne_dir = p.Results.ne_dir;
EEG_stream = p.Results.EEG_stream;
EEG_chan = p.Results.EEG_chan;
EMG_stream = p.Results.EMG_stream;
chan_465 = p.Results.chan_465;
chan_405 = p.Results.chan_405;
chan_ttl_pulse = p.Results.chan_ttl_pulse;
interval = p.Results.interval;
ds_factor_FP = p.Results.ds_factor_FP;
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

if isempty(ne_dir) 
    if isempty(eeg_emg_path)    
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
    if ~strcmp(eeg_emg_path, ne_dir) 
        if isempty(chan_ttl_pulse)
            disp('Please provide TTL pulse channel to sync EEG/EMG with the fp data')
            return
        end
    end
end    

if isempty(save_path)
    save_path = fullfile(parent_dir, strcat(filename, ".mat"));
end

% define the following optional variables
sleep_scores = [];
ne = [];
ne_frequency = nan;
onset_EEG = 0;

%% 2) loading and plotting EEG and EMG raw data

% if EEG/EMG come from FP data
if isempty(exp_file_path)
    fp_data = TDTbin2mat(eeg_emg_path);
    while ~isfield(fp_data.streams, EEG_stream)
        disp(['Invalid EEG stream in the TDT file. You entered', EEG_stream])
        stream_names = fieldnames(fp_data.streams);
        stream_names_char = sprintf('%s, ', stream_names{:});
        message = sprintf('Please type a EEG stream name from the list: {%s}\n', stream_names_char(1:end-2));
        EEG_stream = input(message, 's');
    end

    while ~isfield(fp_data.streams, EMG_stream)
        disp(['Invalid EMG stream in the TDT file. You entered ', EMG_stream])
        stream_names = fieldnames(fp_data.streams);
        stream_names_char = sprintf('%s, ', stream_names{:});
        message = sprintf('Please type a EMG stream name from the list: {%s}\n', stream_names_char(1:end-2));
        EMG_stream = input(message, 's');

    end

    eeg_frequency = fp_data.streams.(EEG_stream).fs; %sampling frequency for EEG signal 
    eeg_data = fp_data.streams.(EEG_stream).data; %EEG signal
    eeg = eeg_data(EEG_chan,:); %add channel (1 or 2)
    emg = fp_data.streams.(EMG_stream).data; %EMG
    total_duration = floor(length(eeg) / eeg_frequency);
    n_seg = ceil(total_duration / 3600 / 12);
    extra_seconds = ceil(n_seg *  3600 * 12 - total_duration);
    %disp(['remainder: ' num2str(remainder)])
    %disp(['n_seg: ' num2str(n_seg)])
    duration_array = 3600 * 12 * ones(1, n_seg); % break into 12-hour segments if necessary
    duration_array(end) = duration_array(end) - extra_seconds;
    bin_filenames = "bin_" + string(1:n_seg);
    video_names = strings(1, n_seg);
    video_dirs = strings(1,n_seg);
else
    Info=loadEXP(exp_file_path,'no');
    bin_filenames = {Info.BinFiles.FileName};
    video_names = {Info.VideosFiles.Files.FileName};
    video_dirs = {Info.VideosFiles.Files.Dir};
    TimeReldebSec = 0; %start extract data from the beginning (first bin)
    % use inf as the end time instead of summing the duration of all bins
    % because duration includes gaps between nins and thus not accurate
    TimeRelEndSec = Inf; 
    [eeg_emg_data, time] = ExtractContinuousData([],Info,[],TimeReldebSec, TimeRelEndSec,[],1);
    emg = eeg_emg_data(1,1:end);
    eeg = eeg_emg_data(2,1:end);
    eeg_frequency = Info.Fs;
    eeg_start_time = [Info.BinFiles.TStart];
    total_duration = length(eeg) / eeg_frequency;
    duration_array = diff(eeg_start_time) * 24 * 3600;
    duration_array = [duration_array total_duration - sum(duration_array)];

end

eeg = single(eeg);
emg = single(emg);
time_eeg = (0:length(eeg)-1)/eeg_frequency;

%% 3) Load FP (fiber_photometry) data (batch II)

if ~isempty(ne_dir)
    if ~strcmp(ne_dir, eeg_emg_path) % if NE data come from different sources
        fp_data = TDTbin2mat(ne_dir); % data is a struct  
    end

    while ~isfield(fp_data.streams, chan_465)
        disp(['Invalid chan_465 in the TDT file. You entered ', chan_465])
        chan_names = fieldnames(fp_data.streams);
        chan_names_char = sprintf('%s, ', chan_names{:});
        message = sprintf('Please type a chan_465 name from the list: {%s}\n', chan_names_char(1:end-2));
        chan_465 = input(message, 's');
    end

    while ~isfield(fp_data.streams, chan_405)
        disp(['Invalid chan_405 in the TDT file. You entered ', chan_405])
        chan_names = fieldnames(fp_data.streams);
        chan_names_char = sprintf('%s, ', chan_names{:});
        message = sprintf('Please type a chan_405 name from the list: {%s}\n', chan_names_char(1:end-2));
        chan_405 = input(message, 's');
    end

    ne_frequency = fp_data.streams.(chan_465).fs; % sampling frequency for NE, one number
    signal_465 = fp_data.streams.(chan_465).data; % hSyn-NE, array 1-D
    signal_405 = fp_data.streams.(chan_405).data; % autofluorescence, array, 1-D
    
    % if different soruces we need to sync it with fp recording
    if ~strcmp(ne_dir, eeg_emg_path)
        % Need to remove FP trace prior to first TTL pulse
        while ~isfield(fp_data.epocs, chan_ttl_pulse)
            disp(['Invalid chan_ttl_pulse in the TDT file. You entered ', chan_ttl_pulse])
            pulse_names = fieldnames(fp_data.epocs);
            pulse_names_char = sprintf('%s, ', pulse_names{:});
            message = sprintf('Please type a chan_ttl_pulse name from the list: {%s}\n', pulse_names_char(1:end-2));
            chan_ttl_pulse = input(message, 's');
        end

        TTL_FP = fp_data.epocs.(chan_ttl_pulse).onset; % TTL_FP is the timestamps
        TTL_gap = diff(TTL_FP) > 5 + 1; % the interval of the pulse is 5 seconds 
        if isempty(find(TTL_gap == 1, 1))
            TTL_onset = TTL_FP(1);  % when TTL pulse train is only started once
        else 
            TTL_onset = TTL_FP(find(TTL_gap==1)+1); % when TTL pulse train is started more than once
        end
        
        first_TTL = TTL_onset(1)*ne_frequency; %sampling point # to start with
        onset_FP = round(first_TTL);
        % remove FP trace prior to first TTL pulse
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
    %[p,~,mu] = polyfit(signal_405(round(mouse{5}*signal_fs)); % for scaling and centering (matlab 2020b and later)
    %controlFit = polyval(p,signal_405,[],mu); % for scaling and centering
    controlFit =  filtfilt(MeanFilter,1,double(controlFit));
    normDat = (signal_465 - controlFit)./controlFit;
    delta_465 = normDat * 100;
    
    % smoothing traces
    ne = filtfilt(MeanFilter,1,double(delta_465));
    
    % downsample NE
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

%% 5) Alingment of EEG recording with respect to FP recording
% if EEG/EMG and fp data come from different sources, align them using the
% TTL pulse.
if ~isempty(ne_dir) && ~strcmp(eeg_emg_path, ne_dir)
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
            % this assumes the "full" setup option for ttl pulse, meaning
            % the ttl pulse was active throughout the fp recording
            onset_EEG = TTL_pulse_time(find(TTL_pulse_time_gap, 1)+1); 
        end    
        onset_EEG_ind = round(onset_EEG*eeg_frequency);
    end
    
    %Cutting EEG/EMG traces leading up to first TTL 
    % Removing first seconds of EEG and EMG raw traces to align with FP trace
    emg = emg(onset_EEG_ind:end);
    eeg = eeg(onset_EEG_ind:end);
    total_duration = floor(length(eeg) / eeg_frequency);
    duration_array(1) = duration_array(1) - round(onset_EEG);
    time_eeg = (0:length(eeg)-1)/eeg_frequency;
    
    if ~isempty(sleep_scores)
        % Remove first seconds of EEG score to align with FP trace
        sleep_scores = sleep_scores(round(onset_EEG+1):end);
    end
end

%% 6) interpolate, extrapolate if there are missing values
% Find the indices of NaN and non-NaN values
nan_indices = isnan(eeg);
if ~isempty(nan_indices)
    non_nan_indices = ~isnan(eeg);
    % Perform interpolation to fill NaN values
    eeg(nan_indices) = interp1(find(non_nan_indices), eeg(non_nan_indices), find(nan_indices), 'linear','extrap');
end

nan_indices = isnan(emg);
if ~isempty(nan_indices)
    non_nan_indices = ~isnan(emg);
    emg(nan_indices) = interp1(find(non_nan_indices), emg(non_nan_indices), find(nan_indices), 'linear','extrap');
end

nan_indices = isnan(ne);
if ~isempty(nan_indices)
    non_nan_indices = ~isnan(ne);
    ne(nan_indices) = interp1(find(non_nan_indices), ne(non_nan_indices), find(nan_indices), 'linear','extrap');
end

 %% 6) plot (optioinal) and save extracted data to .mat file
if show_figure
    if ~isempty(ne_dir)
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

%% 7) segment and (if multiple bins from exp file or if longer than 12 hours) save
num_class = 3;
%disp(['total duration: ', num2str(total_duration)])
%disp(['sleep score len: ', num2str(length(sleep_scores))])
if ~isempty(sleep_scores)
    fill_array = NaN(1, max([0 total_duration - length(sleep_scores)]));
    sleep_scores = [sleep_scores fill_array];
end

if ~isempty(ne)
    min_ne_len = ceil(total_duration * ne_frequency);
    ne_fill_array = NaN(1, max([0 min_ne_len - length(ne)]));
    ne = [ne ne_fill_array];
end


[parent_dir, data_name, ext] = fileparts(save_path);
prev_end  = 0;
n_bins = length(duration_array);
if n_bins > 1
    for i = 1:n_bins
        start_time = prev_end;
        time_end = start_time + floor(duration_array(i));
        prev_end = time_end;
        recording.eeg = eeg(floor(start_time*eeg_frequency+1):floor(time_end*eeg_frequency));
        recording.emg = emg(floor(start_time*eeg_frequency+1):floor(time_end*eeg_frequency));
        if ~isempty(ne)
            recording.ne = ne(floor(start_time*ne_frequency+1):floor(time_end*ne_frequency));
        else
            recording.ne = ne;
        end
        
        if ~isempty(sleep_scores)
            recording.sleep_scores = sleep_scores(start_time+1:time_end);
        end
        
        recording.start_time = start_time;
        recording.num_class = num_class;
        recording.eeg_frequency = eeg_frequency;
        recording.ne_frequency = ne_frequency;
        recording.video_name = video_names{i};
        recording.video_path = fullfile(video_dirs{i}, video_names{i});
        if i == 1
            recording.video_start_time = round(onset_EEG);
        else
            recording.video_start_time = 0;
        end
        bin_filename = bin_filenames{i};
        [~, bin_save_name, ~] = fileparts(bin_filename);
        save_path = fullfile(parent_dir, strcat(data_name, '_', bin_save_name, ".mat"));
        save(save_path, "-struct","recording")
    end
else
    start_time = 0;
    video_start_time = round(onset_EEG);
    video_name = video_names{1};
    video_path = fullfile(video_dirs{1}, video_names{1});
    save(save_path, "eeg", "emg", "ne", "sleep_scores", "start_time", "video_start_time", "num_class", "eeg_frequency", "ne_frequency", "video_name", "video_path")
end
end