% Load the EEG data using EEGLAB
close all;
clear all;


frequncey='500'
stimduration='LLR'
str=['*LLR','*500hz','*.set']

sub_wise=1;
freq_wise =1


filedir='E:\data\epoched';
matfiles = dir(fullfile(filedir, str));
nfiles = length(matfiles);

%%
baseline_removal = 1;
epoch_rejection  = 1;
epoch_reject_max_amplitude = 150;  % If epoch amplitude is higher than that, epoch is rejected
epoch_reject_min_amplitude = -150; % If epoch amplitude is lower than that, epoch is rejected
eeg_channel_pos  = 1:6;

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
    

  %%                Print channel names
            fprintf('Channel Names:\n');
            for i = 1:length({EEG_in.chanlocs.labels})
                fprintf('%s\n', EEG_in.chanlocs(i).labels);
              
            end

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
        params.freqs = [5 100]; % frequency range in Hz
    else
        params.tlimits = [-100 100]; % time limits in ms
        params.cycles = [1 0.5]; % number of cycles in each frequency bin
        params.freqs = [5 100]; % frequency range in Hz
    end 

%%
%     figure;
%     title(cell2mat([subject_id ' ' stimdur ' ' FREQ ' ' event_names ' ' 'ERPS']))
    for c =1:length(eeg_channel_pos)
     [ersp_in,itc_in,powbase,times_in,freqs_in]= newtimef(EEG_in.data(c,:,:), EEG_in.pnts, params.tlimits, EEG_in.srate, params.cycles,'freqs', params.freqs,'plotersp','off','plotitc','off');
     [ersp_anti,itc_anti,powbase,times_anti,freqs_anti]= newtimef(EEG_anti.data(c,:,:), EEG_anti.pnts, params.tlimits, EEG_anti.srate, params.cycles,'freqs', params.freqs,'plotersp','off','plotitc','off');
      ersp_diff(:,:,c,sub)=ersp_in-ersp_anti   ;
      itc_diff(:,:,c,sub)=abs(itc_in-itc_anti);
     end

   sub=sub+1;
    

end 

%%
if freq_wise
     for c=eeg_channel_pos
        for s=1:sub-1
         ersp_chan(:,:,s)=ersp_diff(:,:,c,s) ;
         itc_chan(:,:,s)=itc_diff(:,:,c,s) ;
        end
       
        ersp_avg(:,:,c)=mean(ersp_chan,3);
        itc_avg(:,:,c)=mean(itc_chan,3);
    
    
    
        figure;
        subplot(3, 3, [2 3 5 6]);
        imagesc(times_in, freqs_in, mean(ersp_chan,3));
        title([EEG_in.chanlocs(c).labels,' ' ,stimduration,' ',frequncey ,' ', 'ERPS Average']);
        axis xy;
        colorbar;
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        % Add frequency response plot to the left side
        subplot(3, 3, [1 4]);
        plot(mean(mean(ersp_chan,3),2),freqs_in, 'LineWidth', 2);
        xlabel('Amplitude');
        ylabel('Freq(hz)');
        title('Frequency Response');
%         % % Add ERP plot to the bottom
%         subplot(3,3,[ 8 9]);
%         plot(times_in,mean(mean(ersp_chan,3),1));
%         xlabel('Time (ms)');
%         ylabel('Amplitude');
%         title('ERP');
    
        %  now for TIC
         figure;
        subplot(3, 3, [2 3 5 6]);
        imagesc(times_in, freqs_in, mean(itc_chan,3));
        title([EEG_in.chanlocs(c).labels,' ' ,stimduration,' ',frequncey ,' ', 'ITC Average']);
        axis xy;
        colorbar;
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        % Add frequency response plot to the left side
        subplot(3, 3, [1 4]);
        plot(mean(mean(itc_chan,3),2),freqs_in, 'LineWidth', 2);
        xlabel('Amplitude');
        ylabel('Freq(hz)');
        title('Frequency Response');
        % % Add ERP plot to the bottom
        subplot(3,3,[ 8 9]);
        plot(times_in,mean(mean(itc_chan,3),1));
        xlabel('Time (ms)');
        ylabel('Amplitude');
        title('ERP');
    
    
    
    
     
     end
end

%%
if sub_wise
    for s=1:sub-1
        for c=eeg_channel_pos
         ersp_chan(:,:,c)=ersp_diff(:,:,c,s) ;
         itc_chan(:,:,c)=itc_diff(:,:,c,s) ;
        end
       
        ersp_avg(:,:,s)=mean(ersp_chan,3);
        itc_avg(:,:,s)=mean(itc_chan,3);
    
    
    
        figure;
        subplot(3, 3, [2 3 5 6]);
        imagesc(times_in, freqs_in, mean(ersp_chan,3));
        title(['sub:',num2str(s),' ' ,stimduration,' ',frequncey ,' ', 'ERPS Average']);
        axis xy;
        colorbar;
         xlabel('Time (s)');
         ylabel('Frequency (Hz)');

        % Add frequency response plot to the left side
        subplot(3, 3, [1 4]);
        plot(mean(mean(ersp_chan,3),2),freqs_in, 'LineWidth', 2);
        xlabel('Amplitude');
        ylabel('Freq(hz)');
        title('Frequency Response');
        % % Add ERP plot to the bottom
        subplot(3,3,[ 8 9]);
        plot(times_in,mean(mean(ersp_chan,3),1));
        xlabel('Time (ms)');
        ylabel('Amplitude');
        title('ERP');
    
        %  now for TIC
         figure;
        subplot(3, 3, [2 3 5 6]);
        imagesc(times_in, freqs_in, mean(itc_chan,3));
        title(['sub:',num2str(s),' ' ,stimduration,' ',frequncey ,' ', 'ITC Average']);
        axis xy;
        colorbar;
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');

        % Add frequency response plot to the left side
        subplot(3, 3, [1 4]);
        plot(mean(mean(itc_chan,3),2),freqs_in, 'LineWidth', 2);
        xlabel('Amplitude');
        ylabel('Freq(hz)');
        title('Frequency Response');
        % % Add ERP plot to the bottom
        subplot(3,3,[ 8 9]);
        plot(times_in,mean(mean(itc_chan,3),1));
        xlabel('Time (ms)');
        ylabel('Amplitude');
        title('ERP');
    
    
    
    
     
    end
end
%%

%% cooment this code if ruuning on matlab 2017

    fn = 'E:\data\results\'
    figHandles = findall(0,'Type','figure'); 
    for ind = 1:numel(figHandles)
        exportgraphics(figHandles(ind), [fn,stimduration,'_',frequncey,'_','ERSP_ITC_results.pdf'], 'Append', true);
    end

