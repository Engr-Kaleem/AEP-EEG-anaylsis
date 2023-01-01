close all;
clear all;
FREQ = '500';                 % used to load from inside data folder
stimdur = 'MLR';  
str=['*',stimdur,'*',FREQ,'*.set']

filedir = 'E:\data\epoched';
matfiles = dir(fullfile(filedir, str));
nfiles = length(matfiles);

for i=1:2:nfiles
    subject_id = matfiles(i).name(1:5);       % used to make a folder inside data folder
    FREQ = matfiles(i).name(11:13);                 % used to make a folder inside data folder
    stimdur = matfiles(i).name(7:9);          % used to make a folder inside data folder
    filename_anti=[filedir,'\',matfiles(i).name]
    filename_in=[filedir,'\',matfiles(i+1).name]
       
    EEG_anti = pop_loadset(filename_anti);
    EEG_in = pop_loadset(filename_in);
    

end