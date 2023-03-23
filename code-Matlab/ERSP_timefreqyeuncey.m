% Load the EEG data using EEGLAB
close all;
clear all;




%%
baseline_removal = 0;
epoch_rejection  = 1;
epoch_reject_max_amplitude = 150;  % If epoch amplitude is higher than that, epoch is rejected
epoch_reject_min_amplitude = -150; % If epoch amplitude is lower than that, epoch is rejected
eeg_channel_pos  = 1:12;



%% Load dataset


str=['*.set']

filedir='E:\data\epoched';
matfiles = dir(fullfile(filedir, str));
nfiles = length(matfiles);

for i=1:nfiles
    parts = split(matfiles(i).name, '_');

    % Extract the individual parts and assign them to variables
    subject = parts{1}; % "Sub"
    session = parts{2}; % 1
    stimdur = parts{3}; % "LLR"
    FREQ = parts{4}; % 500
    cond=parts{6}(1:end-4)
    subject_id=[subject '_' ' '  session]
    filename=[filedir,'\',matfiles(i).name]

    EEG = pop_loadset(filename);
    
    
%     Print channel names
%     fprintf('Channel Names:\n');
%     for i = 1:length({EEG.chanlocs.labels})
%         fprintf('%s\n', EEG.chanlocs(i).labels);
%     end
    
    
%     fprintf('Sampling Rate: %d Hz\n', EEG.srate);
    
    
%     Get the event type names
    event_names = unique({EEG.event.type});
    
%     Print the event type names
%     fprintf('Event type names:\n');
%     for i = 1:length(event_names)
%         fprintf('%s\n', event_names{i});
%     end
    
%     Select channels of interest by name
%     channels = {'FC5', 'FC6', 'C5', 'C6','CP5', 'CP6'};
%     chan_inds = [];
%     for i = 1:length(channels)
%         chan_ind = find(strcmp({EEG.chanlocs.labels}, channels{i}));
%         if ~isempty(chan_ind)
%             chan_inds = [chan_inds, chan_ind];
%         end
%     end
%     EEG = pop_select(EEG, 'channel', chan_inds);


     
    epoch_start = -0.3 % in seconds
    epoch_end = 0.4; % in seconds
    EEG = pop_epoch(EEG, { cell2mat(event_names)}, [epoch_start, epoch_end]);
    




%     Baseline correction
if baseline_removal==1
    baseline_start = -0.1; % in seconds
    baseline_end = 0.0; % in seconds
    EEG = pop_rmbase(EEG, [baseline_start, baseline_end]);
end
    
    
     %% Epoch rejection
    % Apply rejections only to eeg channels : channels
    if epoch_rejection
       % [EEG ,rejected_tr]     = rejectbadepochs(EEG, chan_inds, epoch_reject_max_amplitude, epoch_reject_min_amplitude);
        [EEG ,rejected_tr]     = rejectbadepochs(EEG, eeg_channel_pos, epoch_reject_max_amplitude, epoch_reject_min_amplitude);
    end
    
    
    params = struct();
    if stimdur == 'LLR'
            
        params.tlimits = [-100 300]; % time limits in ms
        params.cycles = [1 0.5]; % number of cycles in each frequency bin
        params.freqs = [5 50]; % frequency range in Hz
    else
        params.tlimits = [-100 100]; % time limits in ms
        params.cycles = [1 0.5]; % number of cycles in each frequency bin
        params.freqs = [5 300]; % frequency range in Hz
    end 
    figure;
    title(cell2mat([subject_id ' ' stimdur ' ' FREQ ' ' event_names ' ' 'ERPS']))
    [ersp,itc,powbase,times,freqs] = newtimef(EEG.data(1,:,:), EEG.pnts, params.tlimits, EEG.srate, params.cycles,'freqs', params.freqs);
 
 %  exportgraphics(gcf,'ersp.pdf', 'Append', true);
 %    close all;
    






end

