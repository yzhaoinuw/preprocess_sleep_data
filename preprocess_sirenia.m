function [] = preprocess_sirenia(varargin)

%% 1) Define args and params
signalToolboxDir = fullfile(matlabroot, 'toolbox', 'signal', 'signal');
% Add this directory at the beginning of the MATLAB search path
addpath(signalToolboxDir, '-begin');
%[header, record] = edfread(edfFile);
p = inputParser;
isTextScalar = @(x) ischar(x) || (isstring(x) && isscalar(x));

paramNames = {'edf_file', 'eeg_col', 'emg_col', 'chan_465', 'chan_405', ...
    'chan_ttl_pulse', 'ds_factor_FP', 'save_path', 'show_figure', 'interval'};
if isempty(varargin) || (isTextScalar(varargin{1}) && any(strcmpi(char(varargin{1}), paramNames)))
    varargin = [{''}, varargin];
end

addRequired(p, 'fp_dir', isTextScalar);

default_edf_file = '';
default_eeg_col = 'EEGEEG1A_B';
default_emg_col = 'EMGEMG';
default_chan_465 = 'x465A';
default_chan_405 = 'x405A';
default_chan_ttl_pulse = 'PC0_';
default_ds_factor_FP = 100;
default_save_path = '';
default_show_figure = false;
default_interval = [];

addParameter(p, 'edf_file', default_edf_file, isTextScalar);
addParameter(p, 'eeg_col', default_eeg_col, isTextScalar);
addParameter(p, 'emg_col', default_emg_col, isTextScalar);
addParameter(p, 'chan_465', default_chan_465, isTextScalar);
addParameter(p, 'chan_405', default_chan_405, isTextScalar);
addParameter(p, 'chan_ttl_pulse', default_chan_ttl_pulse, isTextScalar);
addParameter(p, 'ds_factor_FP', default_ds_factor_FP, @isnumeric);
addParameter(p, 'save_path', default_save_path, isTextScalar);
addParameter(p, 'show_figure', default_show_figure, @islogical);
addParameter(p, 'interval', default_interval, @isnumeric);

% Parse the inputs.
parse(p, varargin{:})

% Access the variables.
fp_dir = char(p.Results.fp_dir);

edf_file = char(p.Results.edf_file);
eeg_col = char(p.Results.eeg_col);
emg_col = char(p.Results.emg_col);
chan_465 = char(p.Results.chan_465);
chan_405 = char(p.Results.chan_405);
chan_ttl_pulse = char(p.Results.chan_ttl_pulse);
ds_factor_FP = p.Results.ds_factor_FP;
save_path = char(p.Results.save_path);
show_figure = p.Results.show_figure;
interval = p.Results.interval;

%edf_file = 'C:\Users\yzhao\matlab_projects\pinnacle\Isabelle_M1_Fi21-4hPostICH_M2Fi1PreICH\M1_Fi2PostICH-M2Fi1Pre.edf';
%fp_dir = 'C:\Users\yzhao\matlab_projects\pinnacle\Isabelle_M1_Fi21-4hPostICH_M2Fi1PreICH';

%% 2) Resolve Sirenia and TDT input paths
[~, ~, input_ext] = fileparts(fp_dir);
if isempty(edf_file) && strcmpi(input_ext, '.edf')
    edf_file = fp_dir;
    fp_dir = '';
end

if isempty(edf_file) && ~isempty(fp_dir)
    edf_files = dir(fullfile(fp_dir, '*.edf'));
    % Check if there is exactly one file with that extension
    if isscalar(edf_files)
        %edf_file = edf_files(1).name;
        edf_file = fullfile(edf_files(1).folder, edf_files(1).name);
    elseif isempty(edf_files)
        warning('No .edf file was found in the directory. Please provide the path to the edf file using the argument "edf_file"');
        return
    else
        edfNames = {edf_files.name};
        message = 'More than one .edf file found in fp_dir. Please type one from the list:\n';
        for k = 1:numel(edfNames)
            message = sprintf('%s\x2022 %s\n', message, edfNames{k});
        end
        
        message = [message 'Type your choice: '];
        edfName = input(message, 's');
        edf_file = fullfile(edf_files(1).folder, edfName);
    end
end

if isempty(edf_file)
    warning('No .edf file was provided. Please provide an EDF file path or a folder containing one EDF file.');
    return
end

if ~isfile(edf_file)
    error('Sirenia EDF file was not found: %s', edf_file)
end

has_tdt_data = false;
if ~isempty(fp_dir)
    if ~isfolder(fp_dir)
        error('TDT folder was not found: %s', fp_dir)
    end
    has_tdt_data = ~isempty(dir(fullfile(fp_dir, '*.tsq')));
    if ~has_tdt_data
        disp('No TDT .tsq files found in fp_dir. Processing Sirenia EDF only.')
    end
end

[parent_dir, edfName, ~] = fileparts(edf_file);
if isempty(save_path)
    save_path = fullfile(parent_dir, strcat(edfName, '.mat'));
end

start_time = 0;
sleep_scores = [];
ne = [];
ne_frequency = nan;
video_name = '';
video_path = '';
video_start_time = 0;

%% 3) Extract EEG/EMG from Sirenia
edf_info = edfinfo(edf_file);
recording_start_datetime = datetime( ...
    [char(edf_info.StartDate) ' ' char(edf_info.StartTime)], ...
    'InputFormat', 'dd.MM.yy HH.mm.ss');
[data_table,annotations] = edfread(edf_file);
timeCol = seconds(data_table.('Record Time'));
timeSteps = diff(timeCol);
assert(all(timeSteps(1:end-1) == timeSteps(1)), 'Sampling rate of EEG/EMG changed during the recording.');

eegCol = data_table.(eeg_col);
eeg = vertcat(eegCol{:});
epoch1 = eegCol{1};
eeg_frequency = numel(epoch1) / timeSteps(1);

emgCol = data_table.(emg_col);
emg = vertcat(emgCol{:});

if has_tdt_data
    pulseAltInd = [true; annotations.Annotations(2:end) ~= annotations.Annotations(1:end-1)];
    pulses = annotations(pulseAltInd, :);
    pulseRiseInd = find(contains(pulses.Annotations, 'TTL: Rise'));
    %pulseFallInd = find(contains(pulses.Annotations, 'TTL: Fall'));
    pulseRise = seconds(pulses.Onset(pulseRiseInd));
    %pulseFall = seconds(pulses.Onset(pulseFallInd));
    if isempty(pulseRise)
        error('No "TTL: Rise" annotation was found in the Sirenia EDF file.')
    end
    pulseOnsetInd = max(1, round(pulseRise(1) * eeg_frequency));
else
    pulseOnsetInd = 1;
end

eeg = eeg(pulseOnsetInd:end);
emg = emg(pulseOnsetInd:end);
recording_start_datetime = recording_start_datetime + ...
    seconds((pulseOnsetInd - 1) / eeg_frequency);
eeg = single(eeg);
emg = single(emg);

timeEEG = (0:length(eeg)-1)/eeg_frequency;
%% 4) Extract FP data from TDT and sync with Sirenia

if has_tdt_data
    fpData = TDTbin2mat(fp_dir); % data is a struct
    %ttlPulseOnsets = fpData.epocs.(chan_ttl_pulse).onset;
    %ttlPulseOffsets = fpData.epocs.(chan_ttl_pulse).offset;
    %ttlOnsetDiff = diff(ttlPulseOnsets);
    %ttlOffsetDiff = diff(ttlPulseOffsets);
    %ttlDiff = ttlPulseOffsets - ttlPulseOnsets;

    ne_frequency = fpData.streams.(chan_465).fs; % sampling frequency for NE, one number
    signal_465 = fpData.streams.(chan_465).data; % hSyn-NE, array 1-D
    signal_405 = fpData.streams.(chan_405).data; % autofluorescence, array, 1-D
    ttlFP = fpData.epocs.(chan_ttl_pulse).onset; % ttlFP is the timestamps
    ttlOnset = ttlFP(1); % just keep it simple and take the first one for now

    onsetFPInd = ttlOnset(1)*ne_frequency; %sampling point # to start with
    onsetFPInd = max(1, round(onsetFPInd));
    % remove FP trace prior to first TTL pulse
    signal_465 = signal_465(onsetFPInd:end);
    signal_405 = signal_405(onsetFPInd:end);

    % Normalize and plot (batch II)
    MeanFilterOrder = 1000; % for smoothing
    MeanFilter = ones(MeanFilterOrder,1)/MeanFilterOrder;

    reg = polyfit(signal_405, signal_465, 1);
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
    ne_frequency = ne_frequency / ds_factor_FP;
    timeNE = (0:length(ne)-1)/ne_frequency;

    % sync video start time
    video_start_time = ttlOnset;
end

%% 5) Trim to interval if specified (interval is time in seconds)
if ~isempty(interval)
    t_start = interval(1);
    t_end = interval(2);
    video_start_time = video_start_time + t_start;

    eeg_mask = (timeEEG >= t_start) & (timeEEG <= t_end);
    interval_start_index = find(eeg_mask, 1, 'first');
    if ~isempty(interval_start_index)
        recording_start_datetime = recording_start_datetime + ...
            seconds(timeEEG(interval_start_index));
    end
    eeg = eeg(eeg_mask);
    emg = emg(eeg_mask);
    timeEEG = timeEEG(eeg_mask);

    if has_tdt_data
        ne_mask = (timeNE >= t_start) & (timeNE <= t_end);
        ne = ne(ne_mask);
        timeNE = timeNE(ne_mask);
    end
end

%% 6) Plot
if show_figure
    figure
    if has_tdt_data
        tiledlayout(3, 1, 'TileSpacing', 'compact');
    else
        tiledlayout(2, 1, 'TileSpacing', 'compact');
    end
    
    a = nexttile;
        plot(timeEEG, emg);
        xlabel('time (s)'); ylabel('EMG (V)');
    b = nexttile;
        plot(timeEEG, eeg);
        xlabel('time (s)'); ylabel('EEG (V)');
    if has_tdt_data
        c = nexttile;
            plot(timeNE, ne);
            xlabel('time (s)'); ylabel('NE');
        linkaxes([a, b, c], 'x');
    else
        linkaxes([a, b], 'x');
    end
end

%% 7) Save data
num_class = 3;
recording_start_time = format_recording_start_time(recording_start_datetime);
save(save_path, "eeg", "emg", "ne", "sleep_scores", "start_time", "video_start_time", "recording_start_time", "num_class", "eeg_frequency", "ne_frequency", "video_name", "video_path")
end

function timestamp = format_recording_start_time(start_datetime)
start_datetime.Format = 'yyyy-MM-dd''T''HH:mm:ss.SSS';
timestamp = char(start_datetime);
end
