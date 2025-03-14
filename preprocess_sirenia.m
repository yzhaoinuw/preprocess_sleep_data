function [] = preprocess_sirenia(varargin)

signalToolboxDir = fullfile(matlabroot, 'toolbox', 'signal', 'signal');
% Add this directory at the beginning of the MATLAB search path
addpath(signalToolboxDir, '-begin');
%[header, record] = edfread(edfFile);
p = inputParser;

addRequired(p, 'fp_dir', @ischar);

default_edf_file = '';
default_eeg_col = 'EEGEEG1A_B';
default_emg_col = 'EMGEMG';
default_chan_465 = 'x465A';
default_chan_405 = 'x405A';
default_chan_ttl_pulse = 'PC0_';
default_ds_factor_FP = 100;
default_save_path = '';
default_show_figure = false;

addParameter(p, 'edf_file', default_edf_file, @ischar);
addParameter(p, 'eeg_col', default_eeg_col, @ischar);
addParameter(p, 'emg_col', default_emg_col, @ischar);
addParameter(p, 'chan_465', default_chan_465, @ischar);
addParameter(p, 'chan_405', default_chan_405, @ischar);
addParameter(p, 'chan_ttl_pulse', default_chan_ttl_pulse, @ischar);
addParameter(p, 'ds_factor_FP', default_ds_factor_FP, @isnumeric);
addParameter(p, 'save_path', default_save_path, @ischar);
addParameter(p, 'show_figure', default_show_figure, @islogical);

% Parse the inputs.
parse(p, varargin{:})

% Access the variables.
fp_dir = p.Results.fp_dir;

edf_file = p.Results.edf_file;
eeg_col = p.Results.eeg_col;
emg_col = p.Results.emg_col;
chan_465 = p.Results.chan_465;
chan_405 = p.Results.chan_405;
chan_ttl_pulse = p.Results.chan_ttl_pulse;
ds_factor_FP = p.Results.ds_factor_FP;
save_path = p.Results.save_path;
show_figure = p.Results.show_figure;

%edf_file = 'C:\Users\yzhao\matlab_projects\pinnacle\Isabelle_M1_Fi21-4hPostICH_M2Fi1PreICH\M1_Fi2PostICH-M2Fi1Pre.edf';
%fp_dir = 'C:\Users\yzhao\matlab_projects\pinnacle\Isabelle_M1_Fi21-4hPostICH_M2Fi1PreICH';

if isempty(edf_file)
    edf_files = dir(fullfile(fp_dir, ['*' '.edf']));
    % Check if there is exactly one file with that extension
    if numel(edf_files) == 1
        %edf_file = edf_files(1).name;
        edf_file = fullfile(edf_files(1).folder, edf_files(1).name);
    elseif numel(edf_files) == 0
        warning('No .edf file was found in the directory. Please provide the path to the edf file using the argument "edf_file"');
    else
        edfNames = {edf_files.name};
        message = 'More than one .edf file found in fp_dir. Please type one from the list:\n';
        for k = 1:length(edfNames)
            message = sprintf('%s\x2022 %s\n', message, edfNames{k});
        end
        
        message = [message 'Type your choice: '];
        edfName = input(message, 's');
        edf_file = fullfile(edf_files(1).folder, edfName);
    end
end

[parent_dir, edfName, ~] = fileparts(edf_file);
if isempty(save_path)
    save_path = fullfile(parent_dir, strcat(edfName, '.mat'));
end

start_time = 0;
sleep_scores = [];
video_name = '';
video_path = '';

%% 1) extract EEG/EMG from Sirenia and snyc with TDT
[data_table,annotations] = edfread(edf_file);
timeCol = seconds(data_table.('Record Time'));
timeSteps = diff(timeCol);
assert(all(timeSteps(1:end-1) == timeSteps(1)), 'Sampling rate of EEG/EMG changed during the recording.');

eegCol = data_table.(eeg_col);
eeg = vertcat(eegCol{:});
epoch1 = eegCol{1};
eeg_frequency = length(epoch1) / timeSteps(1);

emgCol = data_table.(emg_col);
emg = vertcat(emgCol{:});

pulseAltInd = [true; annotations.Annotations(2:end) ~= annotations.Annotations(1:end-1)];
pulses = annotations(pulseAltInd, :);
pulseRiseInd = find(contains(pulses.Annotations, 'TTL: Rise'));
%pulseFallInd = find(contains(pulses.Annotations, 'TTL: Fall'));
pulseRise = seconds(pulses.Onset(pulseRiseInd));
%pulseFall = seconds(pulses.Onset(pulseFallInd));
pulseOnsetInd = round(pulseRise(1) * eeg_frequency);

eeg = eeg(pulseOnsetInd:end);
emg = emg(pulseOnsetInd:end);
eeg = single(eeg);
emg = single(emg);

timeEEG = (0:length(eeg)-1)/eeg_frequency;
%% 2) extract fp data from TDT and sync with Sirenia

fpData = TDTbin2mat(fp_dir); % data is a struct 
%ttlPulseOnsets = fpData.epocs.(chan_ttl_pulse).onset;
%ttlPulseOffsets = fpData.epocs.(chan_ttl_pulse).offset;
%ttlOnsetDiff = diff(ttlPulseOnsets);
%ttlOffsetDiff = diff(ttlPulseOffsets);
%ttlDiff = ttlPulseOffsets - ttlPulseOnsets;

ne_frequency = fpData.streams.(chan_465).fs; % sampling frequency for NE, one number
signal_465 = fpData.streams.(chan_465).data; % hSyn-NE, array 1-D
signal_405 = fpData.streams.(chan_405).data; % autofluorescence, array, 1-D
ttlFP = fpData.epocs.(chan_ttl_pulse).onset; % TTL_FP is the timestamps
ttlOnset = ttlFP(1); % just keep it simple and take the first one for now

onsetFPInd = ttlOnset(1)*ne_frequency; %sampling point # to start with
onsetFPInd = round(onsetFPInd);
% remove FP trace prior to first TTL pulse
signal_465 = signal_465(onsetFPInd:end);
signal_405 = signal_405(onsetFPInd:end);

% 3b) Normalize and plot (batch II)
MeanFilterOrder = 1000; % for smoothing
MeanFilter = ones(MeanFilterOrder,1)/MeanFilterOrder;

reg = polyfit(signal_405(1:end), signal_465(1:end), 1);
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

%% 3) plot
if show_figure
    figure
    a = subplot(3,1,1);
        plot(timeEEG, emg); 
        xlabel('time (s)');
        ylabel('EMG (V)');
    b = subplot(3,1,2);
        plot(timeEEG, eeg); 
        xlabel('time (s)');
        ylabel('EEG (V)');
    c = subplot(3,1,3);
        plot(timeNE, ne); 
        xlabel('time (s)');
        ylabel('NE');
    linkaxes([a, b, c],'x');
end

%% 4) save data
num_class = 3;
video_start_time = round(ttlOnset);
save(save_path, "eeg", "emg", "ne", "sleep_scores", "start_time", "video_start_time", "num_class", "eeg_frequency", "ne_frequency", "video_name", "video_path")