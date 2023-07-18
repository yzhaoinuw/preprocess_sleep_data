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

%
% desipramine (batch II)
    %hSyn-NE2m (mPFC), DIO-Arch (LC)
    EEG_FP_477_saline = {'20210712_EEGFP_1_477_2_486_sal' '477_sal/477_sal_2021-07-12_09-53-03-198.exp' '477_sal/477_sal_sleepscore_1s.xlsx' 'x465A' 'x405A' 'PtC0' 3150 1 (1000:12000) 14 15};
    EEG_FP_484_saline = {'20210716_EEGFP_484_sal' '484_sal/484_sal_2021-07-16_09-25-48-573.exp' '484_sal/484_sal_sleepscore_1s.xlsx' 'BuBo' 'V_Bo' 'PtC0' 3030 0 (1000:1800) 14 15};
    EEG_FP_486_saline = {'20210712_EEGFP_1_477_2_486_sal' '486_sal/486_sal_2021-07-12_09-53-03-198.exp' '486_sal/486_sal_sleepscore_1s.xlsx' 'x465C' 'x405C' 'PtC0' 3200 0 (1000:12000) 14 15};
    EEG_FP_569_sal = {'20220807_DPH_580_569' '569_sal/569_sal_2022-08-07_10-34-40-287.exp' '569_sal/569_sal_sleepscore.xlsx' 'x465C' 'x405C' 'Pu1_' 1305 4 (1:10800) 10 11};
%     EEG_FP_578_sal = {'20220807_DPH_578_584' '578_sal/578_sal_2022-08-07_06-48-50-071.exp' '578_sal/578_sleepscore.xlsx' 'x465A' 'x405A' 'Pu1_' 1470 4 (1:10800) 10 11};
    EEG_FP_580_sal = {'20220810_DPH_580_569' '580_sal/580_sal_2022-08-10_10-24-28-685.exp' '580_sal/850_sal_sleepcsore.xlsx' 'x465A' 'x405A' 'Pu1_' 1385 1 (1:10800) 10 11};
    EEG_FP_588_sal = {'20220809_DPH_588_592' '588_sal/588_sal_2022-08-09_10-31-22-994.exp' '588_sal/588_sal_sleepscore.xlsx' 'x465A' 'x405A' 'Pu1_' 1365 1 (1:10800) 10 11};
    EEG_FP_592_sal = {'20220806_diphenH_588_592' '592_sal/592_sal_2022-08-06_10-37-48-796.exp' '592_sal/592_sal_sleepscore.xlsx' 'x465C' 'x405C' 'Pu1_' 1303 0 (1:10800) 10 11};
%     EEG_FP_600_sal = {'20220806_diphenH_602_600' '600/600_2022-08-06_06-52-56-321.exp' '600 - Copy/600_sal_sleepscore.xlsx' 'x465C' 'x405C' 'Pu1_' 1358 3 (1:10800) 10 11};
    EEG_FP_602_sal = {'20220809_DPH_600_602' '602_sal/602_sal_2022-08-09_06-54-24-717.exp' '602_sal/602_sal_sleepscore.xlsx' 'x465C' 'x405C' 'Pu1_' 1305 0 (1:10800) 10 11};
   
   
    mouse = EEG_FP_486_saline;

%% 2b) Load FP data (batch II)

data = TDTbin2mat(mouse{1});

signal_fs = data.streams.(mouse{4}).fs;

signal_465 = data.streams.(mouse{4}).data; %hSyn-NE
signal_405 = data.streams.(mouse{5}).data; %autofluorescence

% if isfield(data.epocs.PtC0.onset)
%     TTL_FP = data.epocs.PtC0.onset.onset;
% elseif isfield(data.epocs.PtC0.onset)
% end

% removing FP trace prior to first TTL pulse
TTL_FP = data.epocs.(mouse{6}).onset;
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

%mouse{8} = (1:7500);
%reg = polyfit(signal_405, signal_465, 1);
reg = polyfit(signal_405(round(mouse{9}*signal_fs)), signal_465(round(mouse{9}*signal_fs)), 1);
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
%delta465_filt = detrend(delta465_filt);
%delta560_filt = filtfilt(MeanFilter,1,double(delta_560));

% downsampling traces for plotting
ds_factor_FP = 100; % also used for plotting later (section 9b)
% % % ds_delta465_filt = delta465_filt;
ds_delta465_filt = downsample(delta465_filt, ds_factor_FP);
%ds_delta560_filt = downsample(delta560_filt, ds_factor_FP);


% % % % % % ds_sec_signal=sec_signal;
ds_sec_signal = downsample(sec_signal, ds_factor_FP); % for plotting

ds_sec_signal_str = datestr(seconds(ds_sec_signal),'HH:MM:SS.FFF');
ds_sec_signal_datetime = datetime(ds_sec_signal_str,'InputFormat','HH:mm:ss.SSS');

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
EEG_time = (1:length(EEG_rawtrace))/sampling_freq;


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

%% 5) open EEG scoring

time_correction = mouse{8}; % NB! If there is a systematic time lag between EEG/EMG traces and scoring adjust for it by seconds here
EEG_sleepscore = xlsread(mouse{3});

% Create binary vectors for sleep stages
%Awake
wake_onset = rmmissing(EEG_sleepscore(:, 2)); % onset of each wake bout in sec (NaNs are removed)
wake_duration = rmmissing(EEG_sleepscore(:, 3)); % duration of each wake bout in sec (NaNs are removed)

%Slow-wave sleep
sws_onset = rmmissing(EEG_sleepscore(:, 6)); % onset of each SWS bout in sec (NaNs are removed)
duration_sws = rmmissing(EEG_sleepscore(:, 7)); % duration of each SWS bout in sec (NaNs are removed)

%REM
    REM_onset = rmmissing(EEG_sleepscore(:, mouse{10})); % onset of each REM bout in sec (NaNs are removed)
    REM_duration = rmmissing(EEG_sleepscore(:, mouse{11})); % duration of each REM bout in sec (NaNs are removed)
    
% Most EEG scorings don't start at time 0 - which shifts the timeline of the
% scoring relative to the EEG/EMG traces - this is corrected for below
if min([wake_onset(1), sws_onset(1), REM_onset(1)]) ~= 0
    EEG_scoring_onset = min([wake_onset(1), sws_onset(1), REM_onset(1)]); % determines the number of seconds to be subtracted
    wake_onset = wake_onset - EEG_scoring_onset;
    sws_onset = sws_onset - EEG_scoring_onset;
    REM_onset = REM_onset - EEG_scoring_onset;
end

% NB! all EEG/EMG traces are not aligned properly with sleep score (~4 s delay)
wake_onset = wake_onset+time_correction; 
sws_onset = sws_onset+time_correction; 
REM_onset = REM_onset+time_correction; 

% Assumption: For binary vectors index 1 = time 0-1s, index 2= time 1-2 sec, and so forth
wake_binary_vector = zeros([1, (Info.HypnoFiles.Duration)+6]); % #243 gives error - add 5 to length
%wake_binary_vector = zeros([1, (Info.HypnoFiles.Duration)+1]); % vector of zeros matching the length of recording in seconds (+1 for last time interval)
for i=1:length(wake_onset) % making time vector for EEG scoring (frequency = 1Hz)
    t = wake_onset(i)+1; % +1 to put time 0 as index 1
    d = wake_duration(i)-1; % -1 compensates for adding 1
    wake_binary_vector(t:t+d) = 1;
end

sws_binary_vector = zeros([1, (Info.HypnoFiles.Duration)+6]); % #243 gives error - add 5 to length
%sws_binary_vector = zeros([1, (Info.HypnoFiles.Duration)+1]); % vector of zeros matching the length of recording in seconds  (+1 for last time interval)
for i=1:length(sws_onset) % making time vector for EEG scoring (frequency = 1Hz)
    t = sws_onset(i)+1; 
    d = duration_sws(i)-1;
    sws_binary_vector(t:t+d) = 1;
end

%REM_binary_vector = zeros([1, (Info.HypnoFiles.Duration)+160]); % #307 gives error - add 160 to length
REM_binary_vector = zeros([1, (Info.HypnoFiles.Duration)+6]); % vector of zeros matching the length of recording in seconds (+1 for last time interval)
if ~isnan(REM_onset)
    for i=1:length(REM_onset) % making time vector for EEG scoring (frequency = 1Hz)
        t = REM_onset(i)+1;
        d = REM_duration(i)-1;
        REM_binary_vector(t:t+d) = 1;
    end
end




% Time vector for sleep scoring (1 Hz)
sleepscore_time = 0:length(wake_binary_vector)-1; % Should be same length for wake/sws/REM binary vectors


figure;
h(1) = subplot(2,1,1);
    plot_sleep(EEG_time, EMG_rawtrace, sleepscore_time, wake_binary_vector, sws_binary_vector, REM_binary_vector);
    xlabel('time (s)');
    ylabel('EMG (V)');
h(2) = subplot(2,1,2);
    plot_sleep(EEG_time, EEG_rawtrace, sleepscore_time, wake_binary_vector, sws_binary_vector, REM_binary_vector);
    xlabel('time (s)');
    ylabel('EEG (V)');
linkaxes([h(1),h(2)],'x');

% 2-column vectors with on- and offsets for each state
wake_periods = [wake_onset wake_onset+wake_duration];
sws_periods = [sws_onset sws_onset+duration_sws];
REM_periods = [REM_onset REM_onset+REM_duration];


%% 6) Dividing wake bouts into microarousals (MA) and wake w/o MA

MA_maxdur = 15; % maximum duration of microarrousal
MA_idx = find(wake_duration < MA_maxdur);
MA_onset = wake_onset(MA_idx);
MA_duration = wake_duration(MA_idx);
MA_binary_vector = zeros([1, (Info.HypnoFiles.Duration)+6]); % #243 gives error - add 5 to length
%MA_binary_vector = zeros([1, (Info.HypnoFiles.Duration)+1]);
for i=1:length(MA_onset) % making time vector for EEG scoring (frequency = 1Hz)
    t = MA_onset(i)+1;
    d = MA_duration(i)-1;
    MA_binary_vector(t:t+d) = 1;
end

% remove micrarrousal from wake vectors
wake_woMA_onset = wake_onset;
wake_woMA_onset(MA_idx) = [];
wake_woMA_duration = wake_duration;
wake_woMA_duration(MA_idx) = [];
wake_woMA_binary_vector = zeros([1, (Info.HypnoFiles.Duration)+6]); % #243 gives error - add 5 to length
%wake_woMA_binary_vector = zeros([1, (Info.HypnoFiles.Duration)+1]);
for i=1:length(wake_woMA_onset) % making time vector for EEG scoring (frequency = 1Hz)
    t = wake_woMA_onset(i)+1;
    d = wake_woMA_duration(i)-1;
    wake_woMA_binary_vector(t:t+d) = 1;
end

% 2-column vectors with on- and offsets for each state
MA_periods = [MA_onset MA_onset+MA_duration];
wake_woMA_periods = [wake_woMA_onset wake_woMA_onset+wake_woMA_duration];

%{
% figure including MA (yellow)
figure
h(1) = subplot(2,1,1);
    plot_sleep(EEG_time, EMG_rawtrace, sleepscore_time, wake_woMA_binary_vector, sws_binary_vector, REM_binary_vector, MA_binary_vector);
    xlabel('time (s)');
    ylabel('EMG (V)');
h(2) = subplot(2,1,2);
    plot_sleep(EEG_time, EEG_rawtrace, sleepscore_time, wake_woMA_binary_vector, sws_binary_vector, REM_binary_vector, MA_binary_vector);
    xlabel('time (s)');
    ylabel('EEG (V)');
linkaxes([h(1),h(2)],'x');
%}



%% 7) State transitions (uncut vectors)

% Creating one vector with different behaviors represented by unique
% numbers (1=wake, 4=sws, 9=REM, 15=MA) at frequency 1Hz
boutscore_vector = zeros([1, (Info.HypnoFiles.Duration)+6]);

% Here using the unaligned "uncut" vectors
for i=1:length(wake_woMA_onset)
    t = wake_woMA_onset(i)+1;
    d = wake_woMA_duration(i)-1;
    boutscore_vector(t:t+d) = 1; % wake=1
end

for i=1:length(sws_onset)
    t = sws_onset(i)+1;
    d = duration_sws(i)-1;
    boutscore_vector(t:t+d) = 4; % sws=4
end

if ~isnan(REM_onset)
    for i=1:length(REM_onset)
        t = REM_onset(i)+1;
        d = REM_duration(i)-1;
        boutscore_vector(t:t+d) = 9; %REM=9
    end
end

for i=1:length(MA_onset)
    t = MA_onset(i)+1;
    d = MA_duration(i)-1;
    boutscore_vector(t:t+d) = 15; %MA=15
end


% Vectors indicate time of transitions in seconds
transition_sws_wake =  find(diff(boutscore_vector)== -3);
transition_wake_sws =  find(diff(boutscore_vector)== 3);
transition_REM_wake =  find(diff(boutscore_vector)== -8);
transition_sws_MA =  find(diff(boutscore_vector)== 11);
transition_REM_sws =  find(diff(boutscore_vector)== -5);
transition_sws_REM =  find(diff(boutscore_vector)== 5);
transition_REM_MA =  find(diff(boutscore_vector)== 6);



%% 9) Alingment of EEG recording and FP recording

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
EEG_time_cut = (1:length(EEG_rawtrace_cut))/sampling_freq;

% Remove first seconds of EEG score to align with FP trace
wake_binary_vector_cut = wake_binary_vector(round(TTL_EEG_onset+1):end);
sws_binary_vector_cut = sws_binary_vector(round(TTL_EEG_onset+1):end);
REM_binary_vector_cut = REM_binary_vector(round(TTL_EEG_onset+1):end);
MA_binary_vector_cut = MA_binary_vector(round(TTL_EEG_onset+1):end);
wake_woMA_binary_vector_cut = wake_woMA_binary_vector(round(TTL_EEG_onset+1):end);

% Align onset, offset, and duration vectors based on TTL
[wake_onset_cut, wake_offset_cut] = binary_to_OnOff(wake_binary_vector_cut);
wake_duration_cut = wake_offset_cut - wake_onset_cut;

[sws_onset_cut, sws_offset_cut] = binary_to_OnOff(sws_binary_vector_cut);
sws_duration_cut = sws_offset_cut - sws_onset_cut;

if ~isnan(REM_onset)
    [REM_onset_cut, REM_offset_cut] = binary_to_OnOff(REM_binary_vector_cut);
    REM_duration_cut = REM_offset_cut - REM_onset_cut;
else
    REM_onset_cut = NaN;    % in case of no REM bouts
    REM_offset_cut = NaN;
    REM_duration_cut = NaN;
end

    [MA_onset_cut, MA_offset_cut] = binary_to_OnOff(MA_binary_vector_cut);
    MA_duration_cut = MA_offset_cut - MA_onset_cut;

    [wake_woMA_onset_cut, wake_woMA_offset_cut] = binary_to_OnOff(wake_woMA_binary_vector_cut);
    wake_woMA_duration_cut = wake_woMA_offset_cut - wake_woMA_onset_cut;

% Align period arrays according to TTL
wake_periods_cut = [wake_onset_cut wake_offset_cut];
sws_periods_cut = [sws_onset_cut sws_offset_cut];
REM_periods_cut = [REM_onset_cut REM_offset_cut];
    MA_periods_cut = [MA_onset_cut MA_offset_cut];
    wake_woMA_periods_cut = [wake_woMA_onset_cut wake_woMA_offset_cut];

% Alignement of transitions vector:
transition_sws_wake_cut = transition_sws_wake-TTL_EEG_onset;
transition_sws_MA_cut = transition_sws_MA-TTL_EEG_onset;


%% 10b) Plotting all traces and scorings together (batch II)

% plot EEG sleep scoring bouts together with FP data

% Time vector for sleep scoring (1 Hz)
% NB! Here using the cut vector to match the length of binary vectors after being cut to align with FP
% % % % % neeg=round(length(EEG_time_cut)/512);
% % % % % ne=round(length(delta465_filt)/signal_fs);
% % % % % if ne<neeg
% % % % %     T=ne*512;
% % % % %     ds_EEG_rawtrace=EEG_rawtrace_cut(1:T);
% % % % %     ds_EMG_rawtrace=EMG_rawtrace_cut(1:T);
% % % % %     fs_signal = 1:1:length(ds_EEG_rawtrace);
% % % % %     ds_EEG_time = fs_signal/512;
% % % % % elseif neeg<ne
% % % % %     T=round(neeg*signal_fs);
% % % % %     ds_delta465_filt=ds_delta465_filt(1:neeg);
% % % % %     fs_signal = 1:1:length(ds_delta465_filt);
% % % % %     ds_sec_signal = fs_signal/signal_fs;
% % % % % else
% % % % %     display('equal')
% % % % % 
% % % % % end


sleepscore_time_cut = 0:length(wake_binary_vector_cut)-1; % should be same length for wake/sws/REM

fig = figure; % un/comment to include/exclude microarousals (MA)
a = subplot(3,1,1);
%     plot_sleep(ds_sec_signal(1000:end), ds_delta2_filt(1000:end), sleepscore_time_cut, wake_binary_vector_cut, sws_binary_vector_cut, REM_binary_vector_cut);
    plot_sleep(ds_sec_signal, ds_delta465_filt, sleepscore_time_cut, wake_woMA_binary_vector_cut, sws_binary_vector_cut, REM_binary_vector_cut, MA_binary_vector_cut);
    title('NE2m');
b = subplot (3,1,2);
    ds_EEG_time = downsample(EEG_time_cut, 10);
    ds_EMG_rawtrace = downsample(EMG_rawtrace_cut, 10);
    %plot_sleep(EEG_time_cut, EMG_rawtrace_cut, sleepscore_time_cut, wake_binary_vector_cut, sws_binary_vector_cut, REM_binary_vector_cut);
    plot_sleep(ds_EEG_time, ds_EMG_rawtrace, sleepscore_time_cut, wake_woMA_binary_vector_cut, sws_binary_vector_cut, REM_binary_vector_cut, MA_binary_vector_cut);
    xlabel('time (s)');
    ylabel('EMG (V)');
c = subplot(3,1,3);
    ds_EEG_rawtrace = downsample(EEG_rawtrace_cut, 10);
    %plot_sleep(EEG_time_cut, EEG_rawtrace_cut, sleepscore_time_cut, wake_binary_vector_cut, sws_binary_vector_cut, REM_binary_vector_cut);
    plot_sleep(ds_EEG_time, ds_EEG_rawtrace, sleepscore_time_cut, wake_woMA_binary_vector_cut, sws_binary_vector_cut, REM_binary_vector_cut, MA_binary_vector_cut);
    xlabel('time (s)');
    ylabel('EEG (V)');
linkaxes([a,b,c],'x');

h = datacursormode(fig);
    h.UpdateFcn = @DataCursor_custom;
    h.SnapToDataVertex = 'on';
    datacursormode on

 %%
nrows_eeg = fix(numel(EEG_rawtrace_cut)/512); % divide by sampling rate to get the duration in seconds
nrows_emg = fix(numel(EMG_rawtrace_cut)/512);
nrows_ne = fix(numel(delta465_filt)/1017);
nrows = min([nrows_eeg nrows_emg nrows_ne]); % take the lowest duration

trial_eeg = reshape(EEG_rawtrace_cut(1:nrows*512), [nrows, 512]);
trial_emg = reshape(EMG_rawtrace_cut(1:nrows*512), [nrows, 512]);
trial_ne = reshape(delta465_filt(1:nrows*1017), [nrows, 1017]);

save data.mat trial_eeg trial_emg trial_ne

%{
trial_eeg=[];
trial_emg=[];
trial_ne=[];
% specify the duration of time like 1000 sec and put T=1000 
% note: ignore the last 20 seconds to avoid length variation thus the 
%final T would be 1000-20

display('specify T value which is the duration of signal')
T=1000;
for i=1:T-1

    while size(EEG_rawtrace_cut,2)>0
        window1=512;
        window2=1017;
        temp_eeg=EEG_rawtrace_cut(1:window1);
        temp_emg=EMG_rawtrace_cut(1:window1);
        temp_ne=delta465_filt(1:window2);
        trial_eeg=[trial_eeg;temp_eeg];
        trial_emg=[trial_emg;temp_emg];
        trial_ne=[trial_ne;temp_ne];
        EEG_rawtrace_cut(1:window)=[];
        EMG_rawtrace_cut(1:window)=[];
        delta465_filt(1:window2)=[];
    end
end


save data2.mat trial_eeg trial_emg trial_ne
%}
