% Load the EEG data using EEGLAB
close all;
clear all;


frequncey='500'
stimduration='LLR'
str=['*LLR','*500hz','*.set']




filedir='E:\data\epoched';
matfiles = dir(fullfile(filedir, str));
nfiles = length(matfiles);

%%
baseline_removal = 1;
epoch_rejection  = 1;
epoch_reject_max_amplitude = 150;  % If epoch amplitude is higher than that, epoch is rejected
epoch_reject_min_amplitude = -150; % If epoch amplitude is lower than that, epoch is rejected
eeg_channel_pos  = 1:12;

sub=1;
 
for i=1:2:nfiles
    parts = split(matfiles(i).name, '_');

    % Extract the individual parts and assign them to variables
    subject = parts{1}; % "Sub"
    session = parts{2}; % 1
    stimdur = parts{3}; % "LLR"
    FREQ = parts{4}; % 500
    cond=parts{6}(1:end-4)
    subject_id=[subject '_' ' '  session]
    filename1=[filedir,'\',matfiles(i).name]
    filename2=[filedir,'\',matfiles(i+1).name]

    EEG_anti = pop_loadset(filename1);
    EEG_in = pop_loadset(filename2);
    
%%            
        %     Print channel names
        %     fprintf('Channel Names:\n');
        %     for i = 1:length({EEG.chanlocs.labels})
        %         fprintf('%s\n', EEG.chanlocs(i).labels);
        %     end
            
%%    
%     fprintf('Sampling Rate: %d Hz\n', EEG.srate);
    
    
%     Get the event type names
    event_in = unique({EEG_in.event.type});
    event_anti = unique({EEG_anti.event.type});
    

 %% 

%     Print the event type names
%     fprintf('Event type names:\n');
%     for i = 1:length(event_names)
%         fprintf('%s\n', event_names{i});
%     end
  %%  
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
%%

     
    epoch_start = -0.3 % in seconds
    epoch_end = 0.4; % in seconds
    EEG_in = pop_epoch(EEG_in, { cell2mat(event_in)}, [epoch_start, epoch_end]);
    EEG_anti = pop_epoch(EEG_anti, { cell2mat(event_anti)}, [epoch_start, epoch_end]);
    


%%

%     Baseline correction
if baseline_removal==1
    baseline_start = -0.1; % in seconds
    baseline_end = 0.0; % in seconds
    EEG_in = pop_rmbase(EEG_in, [baseline_start, baseline_end]);
    EEG_anti = pop_rmbase(EEG_anti, [baseline_start, baseline_end]);
end
 
    
     %% Epoch rejection
    % Apply rejections only to eeg channels : channels
    if epoch_rejection
       % [EEG ,rejected_tr]     = rejectbadepochs(EEG, chan_inds, epoch_reject_max_amplitude, epoch_reject_min_amplitude);
        [EEG_in ,rejected_tr_in]     = rejectbadepochs(EEG_in, eeg_channel_pos, epoch_reject_max_amplitude, epoch_reject_min_amplitude);
        [EEG_anti ,rejected_tr_anti]     = rejectbadepochs(EEG_anti, eeg_channel_pos, epoch_reject_max_amplitude, epoch_reject_min_amplitude);
    end
    
 %%   
    params = struct();
    if stimdur == 'LLR'
            
        params.tlimits = [-100 300]; % time limits in ms
        params.cycles = [1 0.5]; % number of cycles in each frequency bin
        params.freqs = [5 200]; % frequency range in Hz
    else
        params.tlimits = [-100 100]; % time limits in ms
        params.cycles = [1 0.5]; % number of cycles in each frequency bin
        params.freqs = [5 200]; % frequency range in Hz
    end 

%%
%     figure;
%     title(cell2mat([subject_id ' ' stimdur ' ' FREQ ' ' event_names ' ' 'ERPS']))
     [ersp_in,itc_in,powbase,times_in,freqs_in]= newtimef(EEG_in.data(1,:,:), EEG_in.pnts, params.tlimits, EEG_in.srate, params.cycles,'freqs', params.freqs,'plotersp','off','plotitc','off');
     [ersp_anti,itc_anti,powbase,times_anti,freqs_anti]= newtimef(EEG_anti.data(1,:,:), EEG_anti.pnts, params.tlimits, EEG_anti.srate, params.cycles,'freqs', params.freqs,'plotersp','off','plotitc','off');
 
     ersp_diff(:,:,sub)=ersp_anti-ersp_in   ;
     itc_diff(:,:,sub)=abs(itc_anti-itc_in);

   
    figure;
    
%     surf(times_in, freqs_in, ersp_diff);
%     hold on;
    imagesc(times_in, freqs_in,  ersp_diff(:,:,sub));
    axis xy;
    colorbar;
    title('ERSP Difference between inphase and antiphase of ',subject_id);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
   sub=sub+1;
     %  exportgraphics(gcf,'ersp.pdf', 'Append', true);
     %    close all;


end 
%%
figure;
   
%     surf(times_in, freqs_in, ersp_diff);
%     hold on;
    clim=[-1 1];
    imagesc(times_in, freqs_in, mean(ersp_diff,3));
    axis xy;
    colorbar;
    title([stimduration,' ',frequncey ,' ', 'ERPS Average'])
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    figure
    imagesc(times_in, freqs_in, mean(itc_diff,3));
    axis xy;
    colorbar;
    title([stimduration ,' ',frequncey ,' ','ITC average ']);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');

