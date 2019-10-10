% settings seeg analysis
% this script can be updated by the specific user

% add path where you located the repository sEEG_analysis: https://github.com/dvanblooijs/sEEG_analysis
addpath('git_rep/sEEG_analysis')
% add path where you located the repository fieldtrip: https://github.com/fieldtrip/fieldtrip
addpath('git_rep/fieldtrip/')
ft_defaults

% load seeg with SPES from 1 or more patient
cfg.dataPath = '/Fridge/CCEP';
cfg.sub_labels = {'sub-RESP0773'}; % more patients: {'sub-RESPXXXX','sub-RESPXXXX','sub-RESPXXXX'};
cfg.ses_label = 'ses-1';
cfg.task_label = 'task-SPESclin';
cfg.run_label = {'run-021316'}; % more patients: write down the run of each patient: {'run-xxxxxx','run-xxxxxx','run-xxxxxx'};

% settings for SPES 
cfg.epoch_length = 4; % in seconds, -2:2
cfg.epoch_prestim = 2; % seconds before stimulation
cfg.minstim = 5; % number of stimuli needed to be included in further analysis
cfg.dir = 'yes'; % 'yes' or 'no' --> do you want to separate positive and negative stimulations?
cfg.amp = 'no'; % 'yes' or 'no' --> do you want to separate different stimulation currents?

