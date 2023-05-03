% Analysis ERP
close all;
clear all; clc;

eeglab

%% Parameters

dataset_dir  = 'E:\data\epoched';
nsub=2 %number of subjects

FREQ = '750';                 % used to load from inside data folder
stimdur = 'LLR';
STIM_DURATION_MS = 48;     % in ms
EPOCH_TIME_BEFORE_STIM      = -0.020;   % Time in seconds to take before stim to create epoch
EPOCH_TIME_AFTER_STIM       = 0.100;    % Time in seconds to take after stim to create epoch

% Area under the ERP parameters
t_area_min          = EPOCH_TIME_BEFORE_STIM;% for aaread uder the erp curve
t_area_max          = EPOCH_TIME_AFTER_STIM; % for aaread uder the erp curve





f_max   = 200;      % Max frequency for PSD visualisation
freq_bands = {[15, 20], [20, 25], [25, 30] [30, 35], [35, 40], [40, 45]};
freq_bands_names = {'15-20', '20-25', '25-30', '30-35', '35-40', '40-45'};
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

erp_AOC_mat=zeros(nsub,length(eeg_channel_pos));
% freq_engy_mat=zeros(length(eeg_channel_pos),nsub,3*nbands)
freq_engy_mat=[];

f_min   = 2;        % Min frequency for PSD visualisation
f_max   = 200;      % Max frequency for PSD visualisation







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

   
    erp_AOC_mat(sub,1:length(eeg_channel_pos))=erp_energy_anti-erp_energy_in;
    close all;
    sub=sub+1;
end
%%

bpath = 'E:\data\results\'

fn=[bpath,stimdur,'_','ERP_analysis.xlsx']
% Initialize the row and column counters
col_names=[];
for ind=1:length(eeg_channel_pos)
    ind
    col_names{ind} = EEG_in.chanlocs(ind).labels;
end

sheet=FREQ
% Loop over the tables and write them to the Excel file using writetable
data=zeros(nsub+1,length(eeg_channel_pos));
for i = 0:length(eeg_channel_pos)-1
    % Get the current table and its size
   rs=increment_column(i*length(eeg_channel_pos)+1)
   re=increment_column((i+1)*length(eeg_channel_pos))
%    range=[rs,'1',':',re,num2str(nsub+2)]
range=['A1',':','L',num2str(nsub+1)]

    i
 
    data = erp_AOC_mat;
    T = array2table(data, 'VariableNames', col_names)
    writetable(T, fn, 'Sheet', sheet, 'Range', range);

  
end




