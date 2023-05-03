
close all;
clear all; clc;

eeglab

%% Parameters

dataset_dir  = 'E:\data\epoched';
nsub=2 %number of subjects

FREQ = '500';                 % used to load from inside data folder
stimdur = 'MLR';
STIM_DURATION_MS = 48;     % in ms

EPOCH_TIME_BEFORE_STIM      = 0.050;   % Time in seconds to take before stim to create epoch
EPOCH_TIME_AFTER_STIM       = 0.500;    % Time in seconds to take after stim to create epoch

% excel_out = fullfile(dataset_dir, subject_id, BMLD_TYPE, num2str(FREQ), 'freq_analysis.xlsx');
% excel_out = fullfile(dataset_dir,excel_out = fullfile(dataset_dir, subject_id, BMLD_TYPE, num2str(FREQ), 'freq_analysis.xlsx')'freq_analysis.xlsx')

baseline_removal = 0; 
epoch_rejection  = 1;
epoch_reject_max_amplitude = 150;  % If epoch amplitude is higher than that, epoch is rejected
epoch_reject_min_amplitude = -150; % If epoch amplitude is lower than that, epoch is rejected
% eeg_channel_pos  = 3:10;
eeg_channel_pos  = 1:12;


% Frequency Analysis Parameters
t_min   = 0.050;        % Period of time taken for computing the PSD for each epoch : [t_min, t_max];
t_max   = 0.500;      % Period of time taken for computing the PSD for each epoch : [t_min, t_max];
f_min   = 2;        % Min frequency for PSD visualisation
%f_max   = 200;      % Max frequency for PSD visualisation
f_max   = 500;      % Max frequency for PSD visualisation
nfft    = 4096;     % Number of points for the fft
%freq_bands = {[0.5, 4], [4, 8], [8, 15], [16, 31], [32,50]};
freq_bands = {[15, 20], [20, 25], [25, 30] [30, 35], [35, 40], [40, 45], [45, 50]};
freq_bands_names = {'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7'};

% 
% % Plot parameters
% plot_ci             = 1;        % If 1, plot the 95% confidence interval
% ANTIPHASE_COLOR     = 'r';
% INPHASE_COLOR       = 'g';

nbands=7
t_min          = EPOCH_TIME_BEFORE_STIM;  % for psd
t_max          = EPOCH_TIME_AFTER_STIM;  % for PSD




%%
nfft    = 4096;     % Number of points for the fft
baseline_removal = 0;
epoch_rejection  = 1;
epoch_reject_max_amplitude = 150;  % If epoch amplitude is higher than that, epoch is rejected
epoch_reject_min_amplitude = -150; % If epoch amplitude is lower than that, epoch is rejected
eeg_channel_pos  = 1:12;


% freq_engy_mat=zeros(length(eeg_channel_pos),nsub,3*nbands)
freq_engy_mat=[];

f_min   = 2;        % Min frequency for PSD visualisation
f_max   = 200;      % Max frequency for PSD visualisation




% % Plot parameters
% plot_ci             = 0;        % If 1, plot the 95% confidence interval
% ANTIPHASE_COLOR     = 'r';
% INPHASE_COLOR       = 'g';



%% Load dataset


str=['*',stimdur,'*',FREQ,'*.set']

filedir = 'E:\data\epoched';
matfiles = dir(fullfile(filedir, str));
nfiles = length(matfiles);
sub=1
for i=1:2:nfiles
    disp("anylis of subject")
    sub
    subject_id = matfiles(i).name(1:5);       % used to make a folder inside data folder
    FREQ = matfiles(i).name(11:13);                 % used to make a folder inside data folder
    stimdur = matfiles(i).name(7:9);          % used to make a folder inside data folder
    filename_anti=[filedir,'\',matfiles(i).name]
    filename_in=[filedir,'\',matfiles(i+1).name]
       
    EEG_anti = pop_loadset(filename_anti);
    EEG_in = pop_loadset(filename_in);
    
    %% extract epochs
    EEG_anti = pop_epoch(EEG_anti, {  }, [EPOCH_TIME_BEFORE_STIM, EPOCH_TIME_AFTER_STIM], 'epochinfo', 'yes');
    
    EEG_in = pop_epoch(EEG_in, {  }, [EPOCH_TIME_BEFORE_STIM, EPOCH_TIME_AFTER_STIM], 'epochinfo', 'yes');
    
    
    %% Epoch rejection
    % Apply rejections only to eeg channels : channels
    if epoch_rejection
        [EEG_anti,rejected_anti]    = rejectbadepochs(EEG_anti, eeg_channel_pos, epoch_reject_max_amplitude, epoch_reject_min_amplitude);
        [EEG_in ,rejected_in]     = rejectbadepochs(EEG_in, eeg_channel_pos, epoch_reject_max_amplitude, epoch_reject_min_amplitude);
    end
    total_rejected(sub,1)=rejected_anti;
    total_rejected(sub,2)=rejected_in;
    
    %% Baseline Correction
    if baseline_removal ~= 0
        EEG_anti    = pop_rmbase(EEG_anti, [EEG_anti.times(1), 0]);
        EEG_in      = pop_rmbase(EEG_in, [EEG_in.times(1), 0]);
    end
    

   

 
        
  %% Freq anylsis


    for chan_num = eeg_channel_pos
        % Compute mean Power Spectral Density
        [FFT_anti,pxx_db_mean_anti, f_vect]  = freqanalysis_computemeanpsd(EEG_anti, chan_num, nfft, t_min, t_max);
        [FFT_in,pxx_db_mean_in, f_vect]    = freqanalysis_computemeanpsd(EEG_in, chan_num, nfft, t_min, t_max);
        
        % Compute mean power per frequency band
        freq_band_power_anti        = freqanalysis_computebandenergy(pxx_db_mean_anti, f_vect, freq_bands);
        freq_band_power_in          = freqanalysis_computebandenergy(pxx_db_mean_in, f_vect, freq_bands);
        freq_band_diff=freq_band_power_anti-freq_band_power_in;

        freq_engy_mat(chan_num).data(sub,1:nbands)=freq_band_diff;
           
        channel_name = EEG_in.chanlocs(chan_num).labels;
    end

    sub=sub+1;
end


%%



bpath = 'E:\data\results\'

fn=[bpath,stimdur,'_','freq1_analysis.xlsx']
% Initialize the row and column counters
col_names = {'15-20', '20-25', '25-30', '30-35', '35-40', '40-45', '45-50'};

sheet=FREQ
% Loop over the tables and write them to the Excel file using writetable
data=zeros(nsub+1,nbands);
for i = 0:length(eeg_channel_pos)-1
    % Get the current table and its size
   rs=increment_column(i*nbands+1)
   re=increment_column((i+1)*nbands)
   range=[rs,'1',':',re,num2str(nsub+2)]

    i
    data(1,:)=(i+1);
    data(2:end,:) = freq_engy_mat(i+1).data;

    T = array2table(data, 'VariableNames', col_names)
    xlswrite(fn,data,sheet,range)
%     writetable(T, fn, 'Sheet', sheet, 'Range', range);

  
end






