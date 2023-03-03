close all; 
clear all;

%%   add feildtrip toolbox path
addpath C:\Users\kaleem\Documents\MATLAB\fieldtrip-20220707;
ft_defaults;


cfg = [];
cfg.dataset = 'sub-P001_ses-S001_task-Default_run-001_eeg.xdf';
xdfdata = xdf2fieldtrip(cfg.dataset);
 xdfdata.time{1} = xdfdata.time{1} + str2num(xdfdata.hdr.orig.clock_offsets(1).offset{1}.time);
 
%  cfg=[]; 
%  cfg.dataset='sub-PEEL18_ses-S001_task-Default_run-001_eeg.xdf'; 
%  cfg.trialdef.eventtype = 'Markers';
%  cfg.trialdef.prestim=1; 
%  cfg.trialdef.poststim=1; 
%  
%  cfg=ft_definetrial(cfg); 
 