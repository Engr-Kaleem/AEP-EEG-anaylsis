% Analysis ERP
clear all; clc;

eeglab

%% Parameters

dataset_dir  = 'C:\Users\eignatious\OneDrive - Charles Darwin University\Documents\MATLAB\Final_set_Exps\EEG_AEP_Analysis_Final\Frequency_Analysis_29_4_23\Epoched_data';
%dataset_dir  = 'C:\Users\eignatious\OneDrive - Charles Darwin University\Documents\MATLAB\Final_set_Exps\EEG_AEP_Analysis_Final\Analysis_2\Epoched_data';
subject_id = 'Sub_1';       % used to load from inside data folder
FREQ = 1000;                 % used to load from inside data folder
BMLD_TYPE = 'LLR';          % used to load from inside data folder

EPOCH_TIME_BEFORE_STIM      = 0.050;   % Time in seconds to take before stim to create epoch
EPOCH_TIME_AFTER_STIM       = 0.500;    % Time in seconds to take after stim to create epoch


baseline_removal = 0;
epoch_rejection  = 1;
epoch_reject_max_amplitude = 150;  % If epoch amplitude is higher than that, epoch is rejected
epoch_reject_min_amplitude = -150; % If epoch amplitude is lower than that, epoch is rejected
eeg_channel_pos  = 1:12;

% Area under the ERP parameters
t_area_min          = 0.050;
t_area_max          = 0.500;


%N_SESSIONS = 3;
STIM_DURATION_MS = 48;     % in ms

% Plot parameters
plot_ci             = 0;        % If 1, plot the 95% confidence interval
ANTIPHASE_COLOR     = 'r';
INPHASE_COLOR       = 'g';

session             = 1;       % If session == -1, use the merged dataset, if 1 use the 1st session dataset, ...

%% Load dataset
subject_data_dir = fullfile(dataset_dir, subject_id, BMLD_TYPE, num2str(FREQ));

EEG_anti = pop_loadset([subject_id,'_',BMLD_TYPE,'_',num2str(FREQ),'Hz_epoch_anti.set'], subject_data_dir);
EEG_in = pop_loadset([subject_id,'_',BMLD_TYPE,'_',num2str(FREQ),'Hz_epoch_in.set'], subject_data_dir);


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


% %% ERP Images
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




% Analysis ERP
clear all; clc;


