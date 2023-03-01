% Analysis ERP
close all;
clear all; clc;

eeglab

%% Parameters

dataset_dir  = 'E:\data\epoched';
nsub=4 %number of subjects
EPOCH_TIME_BEFORE_STIM      = -0.020;   % Time in seconds to take before stim to create epoch
EPOCH_TIME_AFTER_STIM       = 0.200;    % Time in seconds to take after stim to create epoch
FREQ = '500';                 % used to load from inside data folder
stimdur = 'MLR';  

% Area under the ERP parameters
t_area_min          = EPOCH_TIME_BEFORE_STIM;% for aaread uder the erp curve
t_area_max          = EPOCH_TIME_AFTER_STIM; % for aaread uder the erp curve





f_max   = 200;      % Max frequency for PSD visualisation

t_min          = EPOCH_TIME_BEFORE_STIM;  % for psd
t_max          = EPOCH_TIME_AFTER_STIM;  % for PSD




%%
nfft    = 4096;     % Number of points for the fft
baseline_removal = 0;
epoch_rejection  = 1;

epoch_reject_min_amplitude = -150; % If epoch amplitude is lower than that, epoch is rejected
eeg_channel_pos  = 1:12;



f_min   = 2;        % Min frequency for PSD visualisation
f_max   = 200;      % Max frequency for PSD visualisation




% Plot parameters
plot_ci             = 0;        % If 1, plot the 95% confidence interval
ANTIPHASE_COLOR     = 'r';
INPHASE_COLOR       = 'g';



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
    
    
   
    
    %% Baseline Correction
    if baseline_removal ~= 0
        EEG_anti    = pop_rmbase(EEG_anti, [EEG_anti.times(1), 0]);
        EEG_in      = pop_rmbase(EEG_in, [EEG_in.times(1), 0]);
    end
    
    close all;
  
        
  %%


    for chan_num = eeg_channel_pos
        % Compute mean Power Spectral Density
        [FFT_anti,pxx_db_mean_anti, f_vect]  = freqanalysis_computemeanpsd(EEG_anti, chan_num, nfft, t_min, t_max);
        [FFT_in,pxx_db_mean_in, f_vect]    = freqanalysis_computemeanpsd(EEG_in, chan_num, nfft, t_min, t_max);
        
                
        % Plot the mean PSD
        channel_name = EEG_in.chanlocs(chan_num).labels;
        figure; hold on;
        plot(f_vect, pxx_db_mean_anti, 'color', ANTIPHASE_COLOR);
        plot(f_vect, pxx_db_mean_in, 'color', INPHASE_COLOR);
        axis tight; xlim([f_min, f_max]);
        xlabel('Frequency (Hz)'); ylabel('Gain (dB)');
        title(['Power Spectral Density - Subject ', subject_id,' - channel : ', channel_name, ' - baseline removed : ', num2str(baseline_removal)]);
        legend({'AntiPhase', 'InPhase'});
       
    end

     % save to pdf
     fn = 'E:\data\results\';  %in this example, we'll save to a temp directory.
 
     figHandles = findall(0,'Type','figure'); 

 
%      for ind = 1:numel(figHandles)
%         saveas(figHandles(ind), [fn,subject_id,'_',stimdur,'_',FREQ,'hz_freqAna.pdf'], 'pdf');
%         saveas(gcf,'MyPDFFileName','pdf');
%      end

      for ind = 1:numel(figHandles)
        exportgraphics(figHandles(ind), [fn,subject_id,'_',stimdur,'_',FREQ,'hz_freqAna.pdf'], 'Append', true);
     end
     close all;
     
     
     

    
  
    sub=sub+1;
end


%%
