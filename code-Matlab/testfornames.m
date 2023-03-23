%% load .mat file
close all; 
clear all;





%%    Path of file  you want to use.


filedir = 'E:\data\AEPdata';
matfiles = dir(fullfile(filedir, '*.mat'));
nfiles = length(matfiles);



 


for i=1:nfiles


% Split the filename into parts using the underscore delimiter
parts = split(matfiles(i).name, '_');

% Extract the individual parts and assign them to variables
subject = parts{1}; % "Sub"
session = parts{2}; % 1
type = parts{3}; % "LLR"
frequency = parts{4}(1:end-4); % 500
subject_id=[subject '_'  session]

% Display the extracted parts
disp(['Subject: ', subject])
disp(['Session: ', num2str(session)])
disp(['Type: ', type])
disp(['Frequency: ', num2str(frequency), ' Hz'])
disp(['Subject: ', subject_id])
end