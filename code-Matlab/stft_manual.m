close all;
clear all;

sub_wise=1;
chan_wise =1
frequncey='500'
stimduration='MLR'
str=['*MLR','*500hz','*.set']

only_selected_channels=1;


fs = 2048; % Hz

% Set the frequency range of interest
fmin = 0; % Hz
fmax = 55; % Hz


filedir='E:\data\epoched';
matfiles = dir(fullfile(filedir, str));
nfiles = length(matfiles);

%%
baseline_removal = 1;
epoch_rejection  = 1;
epoch_reject_max_amplitude = 150;  % If epoch amplitude is higher than that, epoch is rejected
epoch_reject_min_amplitude = -150; % If epoch amplitude is lower than that, epoch is rejected
eeg_channel_pos  = 1:3;

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
    

       
%     Get the event type names
    event_in = unique({EEG_in.event.type});
    event_anti = unique({EEG_anti.event.type});
    

  
%%frex       = logspace(log10(10),log10(EEG.srate/5),20);
times2save = -300:25:500;
basetime   = [-300 -100];
timewin    = 300; % in ms

     
    epoch_start = -0.3 % in seconds
    epoch_end = 0.5; % in seconds
    EEG_in = pop_epoch(EEG_in, { cell2mat(event_in)}, [epoch_start, epoch_end]);
    EEG_anti = pop_epoch(EEG_anti, { cell2mat(event_anti)}, [epoch_start, epoch_end]);

    if stimdur == 'LLR'
                  % Set the window length and overlap
        win_length = 0.25 * fs; % in  mili second window
        overlap = round(0.97 * win_length); % 50% overlap
        tmin=0
        tmax=.300;
    else
         % Set the window length and overlap
        win_length = 0.10 * fs; % in  mili second window
        overlap = round(0.99 * win_length); % 50% overlap
        tmin=0
        tmax=.100;
    end 
   
    

    %%    Baseline correction
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
    









% Calculate the STFT for each channel
n_channels = size(EEG_in.data, 1);
% hannels
    for c=eeg_channel_pos

        for ep=1:min(size(EEG_anti.data,3),size(EEG_in.data,3));
                   
                    [si(:,:,ep), fi, ti]=spectrogram(EEG_in.data(c,:,ep), floor(win_length), overlap,1024,fs);
                    [sa(:,:,ep), fa, ta]=spectrogram(EEG_anti.data(c,:,ep), floor(win_length), overlap,1024,fs);

%                 [si(:,:,ep), fi, ti] = stft(EEG_in.data(c,:,ep), fs, 'Window', hamming(win_length), 'OverlapLength', overlap,'FrequencyRange','onesided' );
%                 [sa(:,:,ep), fa, ta] = stft(EEG_anti.data(c,:,ep), fs, 'Window', hamming(win_length), 'OverlapLength', overlap,'FrequencyRange','onesided');
%                 
 
                
                
        end
    si_a =mean(si,3);
    sa_a=mean(sa,3);
    freq_ind=find(fi<=fmax);
    time_ind=find(ti>=0 & ti<=tmax);
    sdiff(:,:,c,sub)=mag2db(abs(si_a(freq_ind,time_ind)-sa_a(freq_ind,time_ind)));
    % Plot the spectrogram for the current channel
%     figure;
%     imagesc(ti, fi, mag2db(abs(sdiff(:,:,c,sub))));
%     axis xy;
%     colorbar;
%     xlabel('Time (s)');
%     ylabel('Frequency (Hz)');
%     title(['sub:',num2str(sub),'Channel ',EEG_in.chanlocs(c).labels,' ' ,stimduration,' ',frequncey ]);
    end
sub=sub+1;
end

%%



if chan_wise
     for c=eeg_channel_pos
        for s=1:sub-1
         s_chan(:,:,s)=sdiff(:,:,c,s) ;
        
        end
       
        s_avg(:,:,c)=mean(s_chan,3);
        
    
    
        figure;
        imagesc(ti(time_ind)*1000, fi(freq_ind), mean(s_chan,3));
        title([EEG_in.chanlocs(c).labels,' ' ,stimduration,' ',frequncey ,' ', 'STFT Average']);
        axis xy;
        colorbar;
        xlabel('Time (ms)');
        ylabel('Frequency (Hz)');
        % Add frequency response plot to the left side
       
    
    
     
     end
end

%%
if sub_wise
    for s=1:sub-1
        for c=eeg_channel_pos
         s_chan(:,:,c)=sdiff(:,:,c,s) ;
        end
       
        s_avg(:,:,s)=mean(s_chan,3);
        figure;
        imagesc(ti(time_ind)*1000, fi(freq_ind), mean(s_chan,3));
        title(['sub:',num2str(s),' ' ,stimduration,' ',frequncey ,' ', 'stft Average']);
        axis xy;
        colorbar;
        xlabel('Time (ms)');
        ylabel('Frequency (Hz)');

    end
end



%%
% s_all=mean(s_suba,3);
% figure;
% imagesc(ti, fi, mag2db(abs(s_all)));
% axis xy;
% colorbar;
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% title([stimduration ,' ',frequncey ,' ','All subject avg ']);
%%


 fn = 'E:\data\results\'
    figHandles = findall(0,'Type','figure'); 
    for ind = 1:numel(figHandles)
        exportgraphics(figHandles(ind), [fn,stimduration,'_',frequncey,'_','STFT_results.pdf'], 'Append', true);
    end


