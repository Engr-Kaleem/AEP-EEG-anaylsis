% Analysis Frequency
clear all; clc;
close all
warning off

eeglab

%addpath('C:\Users\eignatious\OneDrive - Charles Darwin University\Documents\MATLAB\Final_set_Exps\EEG_AEP_Analysis_Final')
%% Parameters
dataset_dir  = 'D:\Google Drive\Upwork\AEPEEGanaylsis\data\Epoched_data';


subject_id = 'Sub_01';       % used to load from inside data folder
FREQ = 1000;                 % used to load from inside data folder
BMLD_TYPE = 'LLR';          % used to load from inside data folder

EPOCH_TIME_BEFORE_STIM      = 0.020;   % Time in seconds to take before stim to create epoch
EPOCH_TIME_AFTER_STIM       = 0.080;    % Time in seconds to take after stim to create epoch

excel_out = fullfile(dataset_dir, subject_id, BMLD_TYPE, num2str(FREQ), 'freq_analysis.xlsx');

baseline_removal = 0; 
epoch_rejection  = 1;
epoch_reject_max_amplitude = 150;  % If epoch amplitude is higher than that, epoch is rejected
epoch_reject_min_amplitude = -150; % If epoch amplitude is lower than that, epoch is rejected
% eeg_channel_pos  = 3:10;
eeg_channel_pos  = 1:12;


% Frequency Analysis Parameters
t_min   = 0.020;        % Period of time taken for computing the PSD for each epoch : [t_min, t_max];
t_max   = 0.080;      % Period of time taken for computing the PSD for each epoch : [t_min, t_max];
f_min   = 2;        % Min frequency for PSD visualisation
%f_max   = 200;      % Max frequency for PSD visualisation
f_max   = 100;      % Max frequency for PSD visualisation
nfft    = 4096;     % Number of points for the fft
%freq_bands = {[0.5, 4], [4, 8], [8, 15], [16, 31], [32,50]};
freq_bands = {[15, 20], [20, 25], [25, 30] [30, 35], [35, 40], [40, 45]};
freq_bands_names = {'15-20', '20-25', '25-30', '30-35', '35-40', '40-45'};

%N_SESSIONS = 3;
STIM_DURATION_MS = 18;     % in ms

% Plot parameters
plot_ci             = 1;        % If 1, plot the 95% confidence interval
ANTIPHASE_COLOR     = 'r';
INPHASE_COLOR       = 'g';


%% Load dataset 
subject_data_dir = fullfile(dataset_dir, subject_id, BMLD_TYPE, num2str(FREQ));

EEG_anti = pop_loadset([subject_id, '_epoch_anti.set'], subject_data_dir);
EEG_in = pop_loadset([subject_id, '_epoch_in.set'], subject_data_dir);


%% extract epochs
EEG_anti = pop_epoch(EEG_anti, {  }, [EPOCH_TIME_BEFORE_STIM, EPOCH_TIME_AFTER_STIM], 'epochinfo', 'yes');

EEG_in = pop_epoch(EEG_in, {  }, [EPOCH_TIME_BEFORE_STIM, EPOCH_TIME_AFTER_STIM], 'epochinfo', 'yes');


%% Epoch rejection
% Apply rejections only to eeg channels : channels 
if epoch_rejection
    EEG_anti    = rejectbadepochs(EEG_anti, eeg_channel_pos, epoch_reject_max_amplitude, epoch_reject_min_amplitude);
    EEG_in      = rejectbadepochs(EEG_in, eeg_channel_pos, epoch_reject_max_amplitude, epoch_reject_min_amplitude);
end

%% Baseline Correction
if baseline_removal ~= 0
    EEG_anti    = pop_rmbase(EEG_anti, [EEG_anti.times(1), 0]);
    EEG_in      = pop_rmbase(EEG_in, [EEG_in.times(1), 0]);
end


%% Mean frequency spectum
% Compute the power spectral density for each trial and compute the mean 
% The PSD is estimated using fft 
k=1;
var{1}='channel_no';
k=k+1;
for i=1:4
    var{k}=strcat('Peak_',num2str(i),'A1_out_mV');
    var{k+1}=strcat('Peak_',num2str(i),'F1_out_Hz');
    var{k+2}=strcat('Peak_',num2str(i),'PH_out_degrees');
    
    var{k+3}=strcat('Peak_',num2str(i),'A1_in_mV');
    var{k+4}=strcat('Peak_',num2str(i),'F1_in_Hz');
    var{k+5}=strcat('Peak_',num2str(i),'PH_in_degrees');
    k=k+6;
end  

peaks=4;
for chan_num = eeg_channel_pos
    % Compute mean Power Spectral Density
    [FFT_anti,pxx_db_mean_anti, f_vect]  = freqanalysis_computemeanpsd(EEG_anti, chan_num, nfft, t_min, t_max);
    [FFT_in,pxx_db_mean_in, f_vect]    = freqanalysis_computemeanpsd(EEG_in, chan_num, nfft, t_min, t_max);
    
    % Compute mean power per frequency band
    freq_band_power_anti        = freqanalysis_computebandenergy(pxx_db_mean_anti, f_vect, freq_bands);
    freq_band_power_in          = freqanalysis_computebandenergy(pxx_db_mean_in, f_vect, freq_bands);
    i=chan_num;
    freq_band_diff=freq_band_power_anti-freq_band_power_in;
    freq_band_power_anti_all(:, i) = freq_band_power_anti;
    freq_band_power_in_all(:, i) = freq_band_power_in;
    freq_band_diff_all(:, i)=freq_band_diff;
    i = i+1;
       
    % Plot the mean PSD
    channel_name = EEG_in.chanlocs(chan_num).labels;
    figure; hold on;
    plot(f_vect, pxx_db_mean_anti, 'color', ANTIPHASE_COLOR);
    plot(f_vect, pxx_db_mean_in, 'color', INPHASE_COLOR);
    axis tight; xlim([f_min, f_max]);
    % Plot frequency band limits
%     ylims = ylim();
%     for freq_band_i = unique(freq_bands);
%         line([freq_band_i, freq_band_i], ylims, 'color', [0.5, 0.5, 0.5], 'linestyle', ':');
%     end
    xlabel('Frequency (Hz)'); ylabel('Gain (dB)');
    title(['Power Spectral Density - Subject ', subject_id,' - channel : ', channel_name, ' - baseline removed : ', num2str(baseline_removal)]);
    legend({'AntiPhase', 'InPhase'});
    
    % Plot mean power per frequency band
    figure; hold on;
    ax(1) = subplot(3,1,1:2);
    x = 1:length(freq_bands);
    h = bar(x, [freq_band_power_anti, freq_band_power_in]);    
    h(1).FaceColor = 'r';
    h(2).FaceColor = 'g';
%     xlabel('Frequency (Hz)');
    ylabel('Gain (dB)');
    ylim([min(min([freq_band_power_anti, freq_band_power_in])) - 2, max(max([freq_band_power_anti, freq_band_power_in])) + 2]);
    set(gca, 'xtick', x); set(gca, 'xticklabel', freq_bands_names);
    title(['Energy per frequency band - Subject ', subject_id,' - channel : ', channel_name, ' - baseline removed : ', num2str(baseline_removal)]);
    legend({'AntiPhase', 'InPhase'});
    ax(2) = subplot(313);
    bar(x, [freq_band_power_anti - freq_band_power_in], 0.5);
    set(gca, 'xtick', x); set(gca, 'xticklabel', freq_bands_names);
    xlabel('Frequency (Hz)'); ylabel('dB');
    title(['Energy per frequency band - AntiPhase - InPhase - Subject ', subject_id,' - channel : ', channel_name, ' - baseline removed : ', num2str(baseline_removal)]);
end

var = {EEG_anti.chanlocs.labels};
var = var(eeg_channel_pos);

var_anti = strcat (var, '_anti');
var_in = strcat (var, '_in');

var = [var_anti, var_in];

freq_band_power_data = [freq_band_power_anti_all, freq_band_power_in_all];

T1 = array2table(freq_band_power_data);
T1.Properties.VariableNames = var;

T1.freq = freq_bands_names';
%T1 = movevars(T1, 'freq', 'Before', 1);

%writetable(T1,excel_out)


