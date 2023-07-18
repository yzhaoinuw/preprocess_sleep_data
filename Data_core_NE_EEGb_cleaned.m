% Data_core_NE_EEG.m: a script for post-processing vital signs that
% indicate sleep state, written by the CTN Copenhagen crew. Updated 18
% January 2023 by Doug Kelley for compatibility with Unix and paths on
% BlueHive. But missing SEV2mat.m, a helper function for TDTbin2mat.m.

%% 1) Define mouse data
clear all 
close all
clc
% data structure:
    % 1) FP raw data
    % 2) EEG raw data
    % 3) EEG sleep score
    % 4) 465 channel name
    % 5) 560 channel name
    % 6) time of injection (s) - back in arena - according to EEG/EMG recording
    % 7) Hypnogram time correction (s)
    % 8) Interval used for signal normalization (only batch II)

% specify experiment files and channels
EEG_FP_477_saline = {'20210712_EEGFP_1_477_2_486_sal' '477_sal/477_sal_2021-07-12_09-53-03-198.exp' 'x465A' 'x405A' 'PtC0' (1000:12000)};
EEG_FP_486_saline = {'20210712_EEGFP_1_477_2_486_sal' '486_sal/486_sal_2021-07-12_09-53-03-198.exp' '486_sal/486_sal_sleepscore_1s.xlsx' 'x465C' 'x405C' 'PtC0' 3200 0 (1000:12000) 14 15};

% {fiber_photometry_(NE)_data  EEG_EMG_data chanA chanB chanC}

mouse = EEG_FP_477_saline;

%% 2b) Load FP (fiber_photometry) data (batch II)

data = TDTbin2mat(mouse{1}); % data is a struct

signal_fs = data.streams.(mouse{3}).fs; % sampling frequency for NE, one number

signal_465 = data.streams.(mouse{3}).data; % hSyn-NE, array 1-D
signal_405 = data.streams.(mouse{4}).data; % autofluorescence, array, 1-D

% removing FP trace prior to first TTL pulse
TTL_FP = data.epocs.(mouse{5}).onset; %
TTL_gap = diff(TTL_FP) > 5 + 1;
if isempty(find(TTL_gap == 1, 1))
    TTL_onset = TTL_FP(1);  % when TTL pulse train is only started once
else 
    TTL_onset = TTL_FP(find(TTL_gap==1)+1); % when TTL pulse train is started more than once
end

first_TTL = TTL_onset(1)*signal_fs; %sampling point # to start with
onset_FP = first_TTL;

signal_465 = signal_465(onset_FP:end);
signal_405 = signal_405(onset_FP:end);


%% 3b) Normalize and plot (batch II)

MeanFilterOrder = 1000; % for smoothing
MeanFilter = ones(MeanFilterOrder,1)/MeanFilterOrder;

fs_signal = 1:1:length(signal_465);
sec_signal = fs_signal/signal_fs;

reg = polyfit(signal_405(round(mouse{6}*signal_fs)), signal_465(round(mouse{6}*signal_fs)), 1);
a = reg(1);
b = reg(2);
controlFit = a.*signal_405 + b;
%[p,~,mu] = polyfit(signal_405(round(mouse{5}*signal_fs));                      % for scaling and centering (matlab 2020b and later)
%controlFit = polyval(p,signal_405,[],mu);                                      % for scaling and centering
controlFit =  filtfilt(MeanFilter,1,double(controlFit));
normDat = (signal_465 - controlFit)./controlFit;
delta_465 = normDat * 100;


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


% smoothing traces
delta465_filt = filtfilt(MeanFilter,1,double(delta_465));

% downsampling traces for plotting
ds_factor_FP = 100; % also used for plotting later (section 9b)
ds_delta465_filt = downsample(delta465_filt, ds_factor_FP);
ds_sec_signal = downsample(sec_signal, ds_factor_FP); % for plotting

% Plot of the three traces above each other (the index 1000:end removes the
% first second of the recoding for nicer plotting)
figure
plot(ds_sec_signal, ds_delta465_filt)
title('NE2m');

%% 4) loading and plotting EEG and EMG raw data

% Make sure the "ExpToolbox" is added to the matlab path and choose the .exp
% file from the folder downloaded from dropbox, e.g. "PROX1_1_EEG" or "PROX1_2_EEG"

% Add functions to path
% addpath(genpath(['Q:Personal_folders/Mie/EEG data from NH/EEG toolbox']));
addpath('EEGtoolbox'); % consider putting this stuff in a common place

% Import EEG raw data to matlab
Info=loadEXP([mouse{2}],'no');

TimeReldebSec=0; %start extract data from the beginning (first bin)
%TimeRelEndSec=inf; %inf to include all data (until last bin)
TimeRelEndSec=Info.BinFiles.Duration; %inf to include all data (including last bin)

[Data,Time]=ExtractContinuousData([],Info,[],TimeReldebSec, TimeRelEndSec,[],1);

EMG_rawtrace = Data(1,1:end);
EEG_rawtrace = Data(2,1:end);

%time vector using sampling frequency
sampling_freq = Info.Fs;
EEG_time = (0:length(EEG_rawtrace)-1)/sampling_freq;


% Plot of EEG and EMG traces

figure
h(1) = subplot(2,1,1);
    plot(EEG_time, EMG_rawtrace); 
    xlabel('time (s)');
    ylabel('EMG (V)');
h(2) = subplot(2,1,2);
    plot(EEG_time, EEG_rawtrace); 
    xlabel('time (s)');
    ylabel('EEG (V)');
linkaxes([h(1),h(2)],'x');


%% 5) Alingment of EEG recording and FP recording

% TTL pulse from FP
TTL_pulse = Data(3,1:end);
onset_EEG = find(diff(TTL_pulse>1*10^-3));
onset_EEG_time = onset_EEG/sampling_freq;
onset_EEG_time_diff = diff(onset_EEG_time);

TTL_gap_EEG = onset_EEG_time_diff > 5;
if isempty(find(TTL_gap_EEG==1, 1))
    onset_EEG = onset_EEG(1);
else 
    onset_EEG = onset_EEG(find(onset_EEG_time_diff>5)+1);
end

TTL_EEG_onset = onset_EEG/sampling_freq;

%Cutting EEG/EMG traces leading up to first TTL 
% Removing first seconds of EEG and EMG raw traces to align with FP trace
EMG_rawtrace_cut = EMG_rawtrace(round(TTL_EEG_onset*sampling_freq):end);
EEG_rawtrace_cut = EEG_rawtrace(round(TTL_EEG_onset*sampling_freq):end);
EEG_time_cut = (0:length(EEG_rawtrace_cut)-1)/sampling_freq;


figure
a = subplot(2,1,1);
    plot(ds_sec_signal, ds_delta465_filt)
    title('NE2m');
b = subplot(2,1,1);
    plot(EEG_time_cut, EMG_rawtrace_cut); 
    xlabel('time (s)');
    ylabel('EMG (V)');
c = subplot(2,1,2);
    plot(EEG_time_cut, EEG_rawtrace_cut); 
    xlabel('time (s)');
    ylabel('EEG (V)');
linkaxes([h(1),h(2)],'x');

 %% 6) reshaping EEG, EMG, and NE

% divide by sampling rate to get the duration in seconds
nrows_eeg = fix(numel(EEG_rawtrace_cut)/512); 
nrows_emg = fix(numel(EMG_rawtrace_cut)/512);
nrows_ne = fix(numel(delta465_filt)/1017);
nrows = min([nrows_eeg nrows_emg nrows_ne]); % take the smallest duration

trial_eeg = transpose(reshape(EEG_rawtrace_cut(1:nrows*512), [512, nrows]));
trial_emg = transpose(reshape(EMG_rawtrace_cut(1:nrows*512), [512, nrows]));
trial_ne = transpose(reshape(delta465_filt(1:nrows*1017), [1017, nrows]));

save data.mat trial_eeg trial_emg trial_ne
