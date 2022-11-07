close all; 
clear all;

%%   add feildtrip toolbox path
addpath D:\fieldtrip-20220104\fieldtrip-20220104;
ft_defaults;
%%

subject_id = 'Sub_1';       % used to make a folder inside data folder
FREQ = 1000;                 % used to make a folder inside data folder
BMLD_TYPE = 'LLR';          % used to make a folder inside data folder

dirpath='D:\Google Drive\Upwork\AEPEEGanaylsis\data\Epoched_data\'
filename= [dirpath,subject_id,'\',BMLD_TYPE,'\',num2str(FREQ),'\',subject_id,'_',BMLD_TYPE,'_',num2str(FREQ),'_epoched.set']

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
cfg.lpfilter        = 'yes';
cfg.lpfreq          = 150;     % low pass cutoff in hertz
data_all = ft_preprocessing(cfg);
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

%% for ploting 
cfg = [];
load easycapM11.mat    % electrode placement layout   make sure the mat file is in same dir as code
cfg.layout = lay;
cfg.interactive = 'yes';
cfg.showoutline = 'yes';
ft_multiplotER(cfg, EEG_in, EEG_anti)



%%  plot the differnce  between  inpahse and anti phase by subtraction

cfg = [];
cfg.operation = 'subtract'; % define operation
cfg.parameter = 'avg';   
difference = ft_math(cfg, EEG_in, EEG_anti); % fine the differncce
%

% now plot the differnce 

cfg = [];
cfg.layout      = lay;
cfg.interactive = 'yes';
cfg.showoutline = 'yes';
ft_multiplotER(cfg, difference);