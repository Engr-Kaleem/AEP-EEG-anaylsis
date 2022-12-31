close all; 
clear all;

%%   add feildtrip toolbox path
addpath C:\Users\kaleem\Documents\MATLAB\fieldtrip-20220707;
ft_defaults;
%%
i=1;
subject_id = 'Sub_1';       % used to make a folder inside data folder
FREQ = 1000;                 % used to make a folder inside data folder
BMLD_TYPE = 'LLR';          % used to make a folder inside data folder

filedir='E:\data\Epoched_data'
matfiles = dir(fullfile(filedir, '*.set'));
nfiles = length(matfiles);

prestim=0.1
poststim=0.2



% for i 1;nfiles
    filename=[filedir,'\',matfiles(i).name]
%%
    cfg = [];
    cfg.dataset =filename; 
    %%
    % check the events  info
    % cfg.trialdef.eventtype = '?';
    % ft_definetrial(cfg);
    
    
    %  define trails
    cfg.trialdef.eventtype      = 'trigger';
    cfg.trialdef.eventvalue     = [-1 1]; % the values of the stimulus trigger for the three conditions
    % 3 = fully incongruent (FIC), 5 = initially congruent (IC), 9 = fully congruent (FC)
    cfg.trialdef.prestim        = 0.1; % in seconds
    cfg.trialdef.poststim       = 0.2; % in seconds
    cfg = ft_definetrial(cfg);
    
    
    
    %% low pass filtering to remove  high frequncey noise 
    % Fitering options and demaan epochs
    cfg.demean          = 'yes';
    cfg.channel={'Cz', 'CPz', 'FCz', 'Pz', 'FC5', 'FC6', 'C5', 'C6', 'CP5', 'CP6', 'T7', 'T8'}
    data_all = ft_preprocessing(cfg);
%%
   cfg = [];
   cfg.viewmode = 'vertical';
   artfct       = ft_databrowser(cfg, data_all)



    %% to see if the triger and stimulli are in sync but do not  do the low pass filter if you wand t to run this code 
    
    % plot(data_all.time{1},data_all.trial{1}(4,:))
    % hold on
    % plot(data_all.time{1},data_all.trial{10}(13,:))
    % plot(data_all.time{1},data_all.trial{1}(15,:))
    
    
    
    
    %% lets seprate the trails for  inphase and anti phase
     % inphase  trails 
    cfg = [];
    cfg.trials = find(data_all.trialinfo==1);
    EEG_in = ft_timelockanalysis(cfg, data_all); % calculate time locked erp
    % antiphase  trails
    cfg = [];
    cfg.trials = find(data_all.trialinfo==-1);
    EEG_anti = ft_timelockanalysis(cfg, data_all); % calculate time locked erp
    
    
    %%
    figure
plot(EEG_anti.time, EEG_anti.avg);
xlabel('time (s)');
ylabel('amplitude (T)');

figure
plot(EEG_in.time, EEG_in.avg);
xlabel('time (s)');
ylabel('amplitude (T/m)');
    
  %%


  for i = 1:12
   
    erp_energy_anti(i) = trapz(EEG_anti.avg(i).^2);
    erp_energy_in(i)   = trapz(EEG_in.avg(i).^2);
  end

x = 1:12;
h = bar(x, [erp_energy_anti, erp_energy_in]);
h(1).FaceColor = 'r';
h(2).FaceColor = 'g';
xlabel('Channel'); ylabel('ERP energy');
set(gca, 'xtick', x); set(gca, 'xticklabel', channel_names(eeg_channel_pos));
legend({'AntiPhase', 'InPhase'});
    
    
    %% for ploting 
%     cfg = [];
%     load easycapM11.mat    % electrode placement layout   make sure the mat file is in same dir as code
%     cfg.layout = lay;
%     cfg.interactive = 'yes';
%     cfg.showoutline = 'yes';
%     ft_multiplotER(cfg, EEG_in, EEG_anti)
%     
    
    
    %%  plot the differnce  between  inpahse and anti phase by subtraction
    
%     cfg = [];
%     cfg.operation = 'subtract'; % define operation
%     cfg.parameter = 'avg';   
%     difference = ft_math(cfg, EEG_in, EEG_anti); % fine the differncce
    %
    
    % now plot the differnce 
    
%     cfg = [];
%     cfg.layout      = lay;
%     cfg.interactive = 'yes';
%     cfg.showoutline = 'yes';
%     ft_multiplotER(cfg, difference);
% end