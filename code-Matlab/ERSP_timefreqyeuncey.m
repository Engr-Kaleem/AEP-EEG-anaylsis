% Load the EEG data using EEGLAB

filename  = 'E:\data\epoched\Sub_2_LLR_500_epoch_in.set';
EEG = pop_loadset(filename);


% Print channel names
fprintf('Channel Names:\n');
for i = 1:length({EEG.chanlocs.labels})
    fprintf('%s\n', EEG.chanlocs(i).labels);
end


fprintf('Sampling Rate: %d Hz\n', EEG.srate);


% Get the event type names
event_names = unique({EEG.event.type});

% Print the event type names
fprintf('Event type names:\n');
for i = 1:length(event_names)
    fprintf('%s\n', event_names{i});
end

% Select channels of interest by name
channels = {'FC5', 'FC6'};
chan_inds = [];
for i = 1:length(channels)
    chan_ind = find(strcmp({EEG.chanlocs.labels}, channels{i}));
    if ~isempty(chan_ind)
        chan_inds = [chan_inds, chan_ind];
    end
end
% EEG = pop_select(EEG, 'channel', chan_inds);

% Epoch the data by event name
event_name = 'InPhase';
epoch_start = -0.2 % in seconds
epoch_end = 0.3; % in seconds
EEG = pop_epoch(EEG, {event_name}, [epoch_start, epoch_end]);

% Baseline correction
baseline_start = -0.2; % in seconds
baseline_end = 0.0; % in seconds
EEG = pop_rmbase(EEG, [baseline_start, baseline_end]);

% Define the frequency bands of interest
freq_bands = [2 :2:70];
times=epoch_start:1/EEG.srate:epoch_end;
% Compute the ERSP using the newtimef function in EEGLAB
ersp = pop_newtimef(EEG,  'freqs', freq_bands, 'nfreqs', size(freq_bands,1),'timesout', times, ...
     'baseline',[baseline_start baseline_end]);

% Plot the ERSP using the plotmatrix function in EEGLAB
figure;
pop_erspimage(ersp,1, [],[],'');
set(gca,'ytick',[1:length(freq_bands)],'yticklabel',{...
    'Delta (2-8Hz)','Theta (8-13Hz)','Alpha (13-20Hz)',...
    'Beta (20-30Hz)','Gamma1 (30-50Hz)','Gamma2 (50-70Hz)'},'fontsize',14);
colorbar;
