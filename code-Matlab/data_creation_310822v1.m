
%% load .mat file
close all; 
clear all;
subject_id = 'Sub_1';       % used to make a folder inside data folder
FREQ = 1000;                 % used to make a folder inside data folder
BMLD_TYPE = 'LLR';          % used to make a folder inside data folder

%%    Path of file  you want to use.

dirpath='D:\Google Drive\Upwork\AEPEEGanaylsis\data\'
filename= [dirpath,subject_id,'_',BMLD_TYPE,'_',num2str(FREQ),'Hz.mat']
load(filename) 
% load ('D:\Google Drive\Upwork\AEPEEGanaylsis\data\Sub_1_LLR_1000Hz.mat');
%%

eeglab;

% basic setup
srate = 19200;
new_srate = 2048;       % downsampled rate
highpass_lowcutoff = 2;   % lower edge of the frequency pass band (Hz)

% filtering;

%% separate 'y' for EEGLAB

time = y(1, :);
stim = y(3, :);
eeg_data = y([6:17, 4,1,3,5,4], :);


%%  To plot some of rawdata  
for i=7:12
x=y(i,19200:19200*5);
plotfr(x,srate);
end;


%% import in EEGLAB

EEG = pop_importdata('dataformat','array','nbchan',0,'data','eeg_data','setname','raw_data','srate',srate,'pnts',0,'xmin',0);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off');
EEG = eeg_checkset( EEG );


%% add channel info
channels_out = {'Cz', 'CPz', 'FCz', 'Pz', 'FC5', 'FC6', 'C5', 'C6', 'CP5', 'CP6', 'T7', 'T8', 'Trigger','tim','stim_L','stim_r','raw_trig'};
channel_loc = struct('labels', channels_out);
EEG.chanlocs = eeg_checkchanlocs(channel_loc);
EEG = eeg_checkset(EEG);


%% resample

EEG = pop_resample( EEG, new_srate);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off');



%% triggers
% adjust triggers for sharp peaks for ease of detection by EEGLAB

trigger = EEG.data(13,:);

[~, locs] = findpeaks(trigger, 'MinPeakHeight', 2e4, 'MinPeakDistance', 0.2*EEG.srate);

if (length(locs) > 500)     % remove any redundant voltage spike in the beginning of recordings
    locs = locs(2:end);
end

if (length(locs) > 500)     % if still more than 500 epochs, take first 500
    locs = locs(1:500);
end

trigger_A = zeros(size(trigger));
trigger_A(locs) = 1;


trigger(trigger > 0) = 0;
trigger = abs(trigger);

[~, locs] = findpeaks(trigger, 'MinPeakHeight', 2e4, 'MinPeakDistance', 0.2*EEG.srate);

if (length(locs) > 500)     % remove any redundant voltage spike in the beginning of recordings
    locs = locs(2:end);
end

if (length(locs) > 500)     % if still more than 500 epochs, take first 500
    locs = locs(1:500);
end

trigger_B = zeros(size(trigger));
trigger_B(locs) = -1;


trigger = trigger_A + trigger_B;
EEG.data(13,:) = trigger;


%% filtering

EEG = pop_eegfiltnew(EEG, 'locutoff', highpass_lowcutoff, 'channels',[1:12] );
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname','filtered_data','gui','off');


%% extract events
EEG = pop_chanevent(EEG, 13,'edge','leading','edgelen',0,'delchan','off');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );

%% save  the data.

% Create output dir if it does not exist
subject_dataset_dir = fullfile('D:\Google Drive\Upwork\AEPEEGanaylsis\data\Epoched_data', subject_id, BMLD_TYPE, num2str(FREQ));
if ~isdir(subject_dataset_dir)
    mkdir(subject_dataset_dir);
end


pop_saveset(EEG, 'filename', [subject_id,'_',BMLD_TYPE,'_',num2str(FREQ),'_epoched'], 'filepath', subject_dataset_dir);







