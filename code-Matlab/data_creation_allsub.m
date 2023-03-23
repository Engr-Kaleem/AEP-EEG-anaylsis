
%% load .mat file
close all; 
clear all;





%%    Path of file  you want to use.


filedir = 'E:\data\AEPdata';
matfiles = dir(fullfile(filedir, '*.mat'));
nfiles = length(matfiles);



 
%% start eeglab
eeglab;
ADJUST_TRIGGERS = 0;       % in ms. values to add in ms to adjust triggers. as triggers were given when trial started. Trial design is |^noise|#tone+noise|noise|. trigger is given at ^, so we need to move trigger at point # as that's where the tone was played.
srate = 19200;

new_srate = 2048;       % downsampled rate

highpass_lowcutoff = 1;   % lower edge of the frequency pass band (Hz)


for i=1:nfiles

    parts = split(matfiles(i).name, '_');

    % Extract the individual parts and assign them to variables
    subject = parts{1}; % "Sub"
    session = parts{2}; % 1
    stimdur = parts{3}; % "LLR"
    FREQ = parts{4}(1:end-4); % 500
    subject_id=[subject '_'  session]

    filename=[filedir,'\',matfiles(i).name]
    
    load(filename) 
    time = y(1, :);
    stim = y(3, :);
    
    eeg_data = y([6:17, 4], :);
%     eeg_data = y([6:17, 4,3,5,4], :);
    
    
    
    %% import in EEGLAB
    
    EEG = pop_importdata('dataformat','array','nbchan',0,'data','eeg_data','setname','raw_data','srate',srate,'pnts',0,'xmin',0);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off');
    EEG = eeg_checkset( EEG );
    
    
    %% add channel info
    channels_out = {'Cz', 'CPz', 'FCz', 'Pz', 'FC5', 'FC6', 'C5', 'C6', 'CP5', 'CP6', 'T7', 'T8', 'Trigger'};
%       channels_out = {'Cz', 'CPz', 'FCz', 'Pz', 'FC5', 'FC6', 'C5', 'C6', 'CP5', 'CP6', 'T7', 'T8', 'Trigger','stim_L','stim_R','Tr'};
    channel_loc = struct('labels', channels_out);
    EEG.chanlocs = eeg_checkchanlocs(channel_loc);
    
    EEG = eeg_checkset(EEG);
    
    
    %% resample
    
    EEG = pop_resample( EEG, new_srate);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off');
    eeglab redraw;
    
    
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
    
    % notch filter
    %EEG = pop_eegfiltnew(EEG, 'locutoff', 45, 'hicutoff', 55, 'channels',[1:12], 'revfilt',1);
    
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname','filtered_data','gui','off');
    eeglab redraw;
    
    
    %% extract events
    EEG = pop_chanevent(EEG, 13,'edge','leading','edgelen',0,'delchan','off');
    EEG = pop_adjustevents(EEG, 'addms',ADJUST_TRIGGERS);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );
    eeglab redraw;
    
    %% save as two separate datasets
    % Create output dir if it does not exist
    subject_dataset_dir = fullfile('E:\data\epoched\');
    if ~isdir(subject_dataset_dir)
        mkdir(subject_dataset_dir);
    end
    
    EEG_anti = pop_selectevent( EEG, 'type',-1,'renametype','AntiPhase','deleteevents','on');
    
    pop_saveset(EEG_anti, 'filename', [subject_id,'_',stimdur,'_',FREQ,'_epoch_anti'], 'filepath', subject_dataset_dir)
    
    EEG_in = pop_selectevent( EEG, 'type',1,'renametype','InPhase','deleteevents','on');
    
    pop_saveset(EEG_in, 'filename', [subject_id,'_',stimdur,'_',FREQ,'_epoch_in'], 'filepath', subject_dataset_dir);
    
    
   
    
end 