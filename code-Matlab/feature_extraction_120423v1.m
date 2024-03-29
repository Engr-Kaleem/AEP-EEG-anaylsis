close all;
clear all;

sub_wise=1;
chan_wise =1
frequncey='500'
stimduration='LLR'
str=['*LLR','*500hz','*.set']

only_selected_channels=1;


fs = 2048; % Hz

% Set the frequency range of interest
fmin = 0; % Hz
fmax = 55; % Hz
nepochs=5

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
    EEG_data = pop_loadset(filename);
   
    

       
%%     Get the event type names
    event_id = unique({EEG_data.event.type});


    if strcmp(cell2mat(event_id) ,'AntiPhase')
       event_label=1;
    else
        event_label=0;
    end
    event_id
    event_label
  %% 
    

  
%%frex       = logspace(log10(10),log10(EEG.srate/5),20);
    times2save = -300:25:500;
    basetime   = [-300 -100];
    timewin    = 300; % in ms

     
    epoch_start = -0.1 % in seconds
    epoch_end = 0.5; % in seconds
    EEG_data = pop_epoch(EEG_data, { cell2mat(event_id)}, [epoch_start, epoch_end]);
    

    if stimdur == 'LLR'
                  % Set the window length and overlap
        win_length = 0.25 * fs; % in  mili second window
        overlap = round(0.97 * win_length); % 50% overlap
        tmin=0.05
        tmax=.250;
        params.tlimits = [-100 300]; % time limits in ms
        params.cycles = [1 0.5]; % number of cycles in each frequency bin
        params.freqs = [5 55]; % frequency range in Hz
    else
         % Set the window length and overlap
        win_length = 0.10 * fs; % in  mili second window
        overlap = round(0.99 * win_length); % 50% overlap
        tmin=0.2
        tmax=.150;
        params.tlimits = [-100 300]; % time limits in ms
        params.cycles = [1 0.5]; % number of cycles in each frequency bin
        params.freqs = [5 55]; % frequency range in Hz
    end 
   
    

    %%    Baseline correction
    if baseline_removal==1
        baseline_start = -0.1; % in seconds
        baseline_end = 0.0; % in seconds
        EEG_data = pop_rmbase(EEG_data, [baseline_start, baseline_end]);
        
    end
 

     %% Epoch rejection
    % Apply rejections only to eeg channels : channels
    if epoch_rejection
       % [EEG ,rejected_tr]     = rejectbadepochs(EEG, chan_inds, epoch_reject_max_amplitude, epoch_reject_min_amplitude);
        [EEG_data ,rejected_tr_in]     = rejectbadepochs(EEG_data, eeg_channel_pos, epoch_reject_max_amplitude, epoch_reject_min_amplitude);
        
    end
    









% % Calculate the STFT for each channel
% n_channels = size(EEG_data.data, 1);
epochs(i)=size(EEG_data.data,3)
% hannels
       for ep=1:size(EEG_data.data,3);
                   for c=eeg_channel_pos;
                    [s_c(:,:,c), fi, ti]=spectrogram(EEG_data.data(c,:,ep), floor(win_length), overlap,1024,fs);
%                     [ersp(:,:,c),itcc,powbase,times,freqs]= newtimef(EEG_data.data(c,:,ep:ep+10), EEG_data.pnts, params.tlimits, EEG_data.srate, params.cycles,'freqs', params.freqs,'plotersp','off','plotitc','off');
%                     itc(:,:,c)=abs(itcc);
                   end
                   freq_ind=find(fi<=fmax);
                   time_ind=find(ti>= tmin & ti<=tmax);
                   s=mean(abs(s_c),3);
                   smat(:,:,ep,i)=mag2db(s(freq_ind,time_ind));
%                    erspmat(:,:,ep,i) =mean(ersp,3)  ;      
%                    itcmat(:,:,ep,i)=mean(abs(itc),3);
                   labels_stft(i,ep)=event_label;
                          
                
       end
       tlabels(i)=floor(epochs(i)/5);
       extra_epochs=mod(epochs(i),5);
       epindx=1;
       for ep=1:nepochs:(size(EEG_data.data,3)-extra_epochs);
                   for c=eeg_channel_pos;
%                    
                    [ersp(:,:,c),itcc,powbase,times,freqs]= newtimef(EEG_data.data(c,:,ep:ep+nepochs-1), EEG_data.pnts, params.tlimits, EEG_data.srate, params.cycles,'freqs', params.freqs,'plotersp','off','plotitc','off');
                   
                    itc(:,:,c)=abs(itcc);
                   end
%                    freq_ind=find(fi<=fmax);
                    
                   time_ind=find(times>=(tmin*1000) & times<=(tmax*1000));
                   erspmat(:,:,epindx,i) =mean(ersp(:,time_ind),3)  ;      
                   itcmat(:,:,epindx,i)=mean(abs(itc(:,time_ind)),3);
                   labels(i,epindx)=event_label;
                   epindx=epindx+1;      
                
        end

end
%%
save('features.mat', 'smat', 'erspmat','itcmat','labels','labels_stft','epochs','tlabels'); % save both variables to a file

