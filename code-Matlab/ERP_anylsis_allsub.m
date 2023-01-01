% Analysis ERP
close all;
clear all; clc;

eeglab

%% Parameters

dataset_dir  = 'E:\data\epoched';
nsub=4
% subject_id = 'Sub_4';       % used to load from inside data folder
FREQ = '500';                 % used to load from inside data folder
stimdur = 'LLR';          % used to load from inside data folder

EPOCH_TIME_BEFORE_STIM      = -0.020;   % Time in seconds to take before stim to create epoch
EPOCH_TIME_AFTER_STIM       = 0.100;    % Time in seconds to take after stim to create epoch



f_max   = 100;      % Max frequency for PSD visualisation
nfft    = 4096;     % Number of points for the fft
%freq_bands = {[0.5, 4], [4, 8], [8, 15], [16, 31], [32,50]};
freq_bands = {[15, 20], [20, 25], [25, 30] [30, 35], [35, 40], [40, 45]};
freq_bands_names = {'15-20', '20-25', '25-30', '30-35', '35-40', '40-45'};


baseline_removal = 0;
epoch_rejection  = 1;
epoch_reject_max_amplitude = 150;  % If epoch amplitude is higher than that, epoch is rejected
epoch_reject_min_amplitude = -150; % If epoch amplitude is lower than that, epoch is rejected
eeg_channel_pos  = 1:12;

% Area under the ERP parameters
t_area_min          = EPOCH_TIME_BEFORE_STIM;
t_area_max          = EPOCH_TIME_AFTER_STIM;

t_min          = EPOCH_TIME_BEFORE_STIM;
t_max          = EPOCH_TIME_AFTER_STIM;
f_min   = 2;        % Min frequency for PSD visualisation
f_max   = 200;      % Max frequency for PSD visualisation


%N_SESSIONS = 3;
STIM_DURATION_MS = 18;     % in ms

% Plot parameters
plot_ci             = 0;        % If 1, plot the 95% confidence interval
ANTIPHASE_COLOR     = 'r';
INPHASE_COLOR       = 'g';

session             = 1;       % If session == -1, use the merged dataset, if 1 use the 1st session dataset, ...

%% Load dataset

REQ = '500';                 % used to load from inside data folder
stimdur = 'MLR';  
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
        EEG_anti    = rejectbadepochs(EEG_anti, eeg_channel_pos, epoch_reject_max_amplitude, epoch_reject_min_amplitude);
        EEG_in      = rejectbadepochs(EEG_in, eeg_channel_pos, epoch_reject_max_amplitude, epoch_reject_min_amplitude);
    end
    
    %% Baseline Correction
    if baseline_removal ~= 0
        EEG_anti    = pop_rmbase(EEG_anti, [EEG_anti.times(1), 0]);
        EEG_in      = pop_rmbase(EEG_in, [EEG_in.times(1), 0]);
    end
    
    %% Plot all trials and ERP
    % for chan_num = 1:EEG_anti.nbchan
    %    plot_alltrialanderp(EEG_anti, chan_num, baseline_removal);
    % end
    %
    % for chan_num = 1:EEG_anti.nbchan
    %    plot_alltrialanderp(EEG_in, chan_num, baseline_removal);
    % end
    close all;
    %% Plot ERPs for AntiPhase and InPhase condition
    for chan_num = eeg_channel_pos
        channel_name = EEG_anti.chanlocs(chan_num).labels;
        figure;
        ax = subplot(1,1,1);
        erp_subplot(EEG_anti, chan_num, plot_ci, ANTIPHASE_COLOR);
        erp_subplot(EEG_in, chan_num, plot_ci, INPHASE_COLOR);
        axis tight;
        plot([0, 0], ylim, 'k'); plot(xlim, [0, 0], 'k--');
        plot([STIM_DURATION_MS, STIM_DURATION_MS], ylim, 'Color', [0.1, 0.1, 0.1]);
        xlabel('Time (ms)'); ylabel('Amplitude (uV)');
        %title(['ERPs - ', 'Subect 2', ' - channel : ', channel_name, ' - baseline removed : ', num2str(baseline_removal)]);
        title(['ERPs - Subject ', subject_id, ' - channel : ', channel_name, ' - baseline removed : ', num2str(baseline_removal)]);
        if plot_ci
            legend({'Antiphase', 'Antiphase 95% CI', 'HomoPhase', 'HomoPhase 95% CI'});
        else
            legend({'Antiphase', 'HomoPhase'});
        end
    end
%%
     

     fn = 'E:\data\results\';  %in this example, we'll save to a temp directory.
 
     figHandles = findall(0,'Type','figure'); 

 
     for ind = 1:numel(figHandles)
        exportgraphics(figHandles(ind), [fn,subject_id,'_',stimdur,'_',FREQ,'hz_ERP.pdf'], 'Append', true);
     end
     close all;
        
  %%


%   tftopo(EEG_anti,EEG_anti.times,allfreqs(:,:,1),'mode','ave','limits', …
% [nan nan nan 35 -1.5 1.5],'signifs', allerspboot, 'sigthresh', [6], 'timefreqs', ...
% [400 8; 350 14; 500 24; 1050 11], 'chanlocs', EEG.chanlocs);

  %% Mean frequency spectum
% Compute the power spectral density for each trial and compute the mean 
% The PSD is estimated using fft 
%     k=1;
%     var{1}='channel_no';
%     k=k+1;
%     for i=1:4
%         var{k}=strcat('Peak_',num2str(i),'A1_out_mV');
%         var{k+1}=strcat('Peak_',num2str(i),'F1_out_Hz');
%         var{k+2}=strcat('Peak_',num2str(i),'PH_out_degrees');
%         
%         var{k+3}=strcat('Peak_',num2str(i),'A1_in_mV');
%         var{k+4}=strcat('Peak_',num2str(i),'F1_in_Hz');
%         var{k+5}=strcat('Peak_',num2str(i),'PH_in_degrees');
%         k=k+6;
%     end  
%     
%     peaks=4;
%     for chan_num = eeg_channel_pos
%         % Compute mean Power Spectral Density
%         [FFT_anti,pxx_db_mean_anti, f_vect]  = freqanalysis_computemeanpsd(EEG_anti, chan_num, nfft, t_min, t_max);
%         [FFT_in,pxx_db_mean_in, f_vect]    = freqanalysis_computemeanpsd(EEG_in, chan_num, nfft, t_min, t_max);
%         
%         % Compute mean power per frequency band
%         freq_band_power_anti        = freqanalysis_computebandenergy(pxx_db_mean_anti, f_vect, freq_bands);
%         freq_band_power_in          = freqanalysis_computebandenergy(pxx_db_mean_in, f_vect, freq_bands);
%         i=chan_num;
%         freq_band_diff=freq_band_power_anti-freq_band_power_in;
%         freq_band_power_anti_all(:, i) = freq_band_power_anti;
%         freq_band_power_in_all(:, i) = freq_band_power_in;
%         freq_band_diff_all(:, i)=freq_band_diff;
%         i = i+1;
%            
%         % Plot the mean PSD
%         channel_name = EEG_in.chanlocs(chan_num).labels;
%         figure; hold on;
%         plot(f_vect, pxx_db_mean_anti, 'color', ANTIPHASE_COLOR);
%         plot(f_vect, pxx_db_mean_in, 'color', INPHASE_COLOR);
%         axis tight; xlim([f_min, f_max]);
%         % Plot frequency band limits
%     %     ylims = ylim();
%     %     for freq_band_i = unique(freq_bands);
%     %         line([freq_band_i, freq_band_i], ylims, 'color', [0.5, 0.5, 0.5], 'linestyle', ':');
%     %     end
%         xlabel('Frequency (Hz)'); ylabel('Gain (dB)');
%         title(['Power Spectral Density - Subject ', subject_id,' - channel : ', channel_name, ' - baseline removed : ', num2str(baseline_removal)]);
%         legend({'AntiPhase', 'InPhase'});
%         
%         % Plot mean power per frequency band
%         figure; hold on;
%         ax(1) = subplot(3,1,1:2);
%         x = 1:length(freq_bands);
%         h = bar(x, [freq_band_power_anti, freq_band_power_in]);    
%         h(1).FaceColor = 'r';
%         h(2).FaceColor = 'g';
%     %     xlabel('Frequency (Hz)');
%         ylabel('Gain (dB)');
%         ylim([min(min([freq_band_power_anti, freq_band_power_in])) - 2, max(max([freq_band_power_anti, freq_band_power_in])) + 2]);
%         set(gca, 'xtick', x); set(gca, 'xticklabel', freq_bands_names);
%         title(['Energy per frequency band - Subject ', subject_id,' - channel : ', channel_name, ' - baseline removed : ', num2str(baseline_removal)]);
%         legend({'AntiPhase', 'InPhase'});
%         ax(2) = subplot(313);
%         bar(x, [freq_band_power_anti - freq_band_power_in], 0.5);
%         set(gca, 'xtick', x); set(gca, 'xticklabel', freq_bands_names);
%         xlabel('Frequency (Hz)'); ylabel('dB');
%         title(['Energy per frequency band - AntiPhase - InPhase - Subject ', subject_id,' - channel : ', channel_name, ' - baseline removed : ', num2str(baseline_removal)]);
%     end
% 
%    
%      fn = 'E:\data\results\';  %in this example, we'll save to a temp directory.
%  
%      figHandles = findall(0,'Type','figure'); 
% 
%  
%      for ind = 1:numel(figHandles)
%         exportgraphics(figHandles(ind), [fn,subject_id,'_',stimdur,'_',FREQ,'hz_freqAna.pdf'], 'Append', true);
%      end
%      close all;
     
     
     
     %% ERP Images
    % for chan_num = eeg_channel_pos
    %     channel_name = EEG_in.chanlocs(chan_num).labels;
    %     
    %         figure;
    %         data_anti = squeeze(EEG_anti.data(chan_num, :, :));
    %         title_str = ['ERP Image - Subject ', subject_id,' AntiPhase - channel : ', channel_name, ' - baseline removed : ', num2str(baseline_removal)];
    %         erpimage(data_anti, [], EEG_anti.times, title_str, 1, 0, 'erp', 3, 'cbar', 'on');
    %     
    %         figure;
    %         data_in = squeeze(EEG_in.data(chan_num, :, :));
    %         title_str = ['ERP Image - Subject ', subject_id,' InPhase - channel : ', channel_name, ' - baseline removed : ', num2str(baseline_removal)];
    %         erpimage(data_in, [], EEG_in.times, title_str, 1, 0, 'erp', 3, 'cbar', 'on');
    % end
    
    
    %% Area under the ERP²
    % Compute the area under the squared ERP for both condition, for all
    % channels.
    % Compute the area between t_area_min and t_area_max
    
    n_eeg_channels  = length(eeg_channel_pos);
    erp_energy_anti = zeros(n_eeg_channels, 1);
    erp_energy_in   = zeros(n_eeg_channels, 1);
    t_sel_ind       = (EEG_in.times > t_area_min*1000) & (EEG_in.times < t_area_max*1000);
    for i = 1:n_eeg_channels
        erp_i_anti = mean(squeeze(EEG_anti.data(eeg_channel_pos(i), :, :)), 2);
        erp_i_in = mean(squeeze(EEG_in.data(eeg_channel_pos(i), :, :)), 2);
        erp_energy_anti(i) = trapz(erp_i_anti(t_sel_ind).^2);
        erp_energy_in(i)   = trapz(erp_i_in(t_sel_ind).^2);
    end
    
    channel_names = {EEG_anti.chanlocs.labels};
    figure; hold on;
    x = 1:n_eeg_channels;
    h = bar(x, [erp_energy_anti, erp_energy_in]);
    h(1).FaceColor = 'r';
    h(2).FaceColor = 'g';
    xlabel('Channel'); ylabel('ERP energy');
    set(gca, 'xtick', x); set(gca, 'xticklabel', channel_names(eeg_channel_pos));
    legend({'AntiPhase', 'InPhase'});
    %title_str = ['ERP mean energy between ', num2str(t_area_min),'s and ', num2str(t_area_max),'s - Subject ', '2', ' - baseline removed : ', num2str(baseline_removal)];
    title_str = ['ERP mean energy between ', num2str(t_area_min),'s and ', num2str(t_area_max),'s - Subject ', subject_id, ' - baseline removed : ', num2str(baseline_removal)];
    title(title_str);
    saveas(gcf,[fn,subject_id,'_',stimdur,'_',FREQ,'hz.png'])
    close all;
    sub=sub+1;
end





