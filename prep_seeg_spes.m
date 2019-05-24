% ER analysis in seeg
% this code is used to load seeg files in BIDS and analyze the early
% responses
% author: D van Blooijs
% date: May 2019

addpath('Desktop/git_rep/sEEG_analysis')
addpath('Desktop/git_rep/fieldtrip/')
ft_defaults

%% load seeg with SPES from 1 patients

dataPath = '/Fridge/CCEP';
sub_labels = { 'RESP0773'};
ses_label = '1';
task_label = 'SPESclin';

i=1;
D = dir(fullfile(dataPath,['sub-' sub_labels{i}],['ses-' ses_label],'ieeg',...
    ['sub-' sub_labels{i} '_ses-' ses_label '_task-' task_label ,'_run-*_ieeg.eeg']));

dataName = fullfile(D(1).folder, D(1).name);

ccep_data = ft_read_data(dataName,'dataformat','brainvision_eeg');
ccep_header = ft_read_header(dataName);

% load events
D = dir(fullfile(dataPath,['sub-' sub_labels{i}],['ses-' ses_label],'ieeg',...
    ['sub-' sub_labels{i} '_ses-' ses_label '_task-' task_label ,'_run-*_events.tsv']));

eventsName = fullfile(D(1).folder, D(1).name);

cc_events = readtable(eventsName,'FileType','text','Delimiter','\t');

% load channels
D = dir(fullfile(dataPath,['sub-' sub_labels{i}],['ses-' ses_label],'ieeg',...
    ['sub-' sub_labels{i} '_ses-' ses_label '_task-' task_label ,'_run-*_channels.tsv']));

channelsName = fullfile(D(1).folder, D(1).name);

cc_channels = readtable(channelsName,'FileType','text','Delimiter','\t');
channels_incl = strcmp(cc_channels.status_description,'included');

data = ccep_data(channels_incl,:);
ch = cc_channels.name(channels_incl);

%% epoch files in 2spre-2spost stimulus 
epoch_length = 4; % in seconds, -2:2
epoch_prestim = 2;
fs = ccep_header.Fs;

stimevent = (strcmp(cc_events.trial_type,'electrical_stimulation'));
stim1mA = double(strcmp(cc_events.electrical_stimulation_current,'0.001'));
stim2mA = double(strcmp(cc_events.electrical_stimulation_current,'0.002'));
stimcur = stim1mA + 2*stim2mA ;
stimchans = cc_events.electrical_stimulation_site(stimevent);
data_epoch = zeros(size(data,1),sum(stimevent(:,1)),round(epoch_length*fs));

tt = (1/fs:1/fs:epoch_length) - epoch_prestim;

for elec = 1:size(data,1) % for all channels
    n=1;
    for ll = 1:sum(stimevent) % for all epochs
        stimnum = find(stimevent ==1);
        data_epoch(elec,n,:) = data(elec,cc_events.sample_start(stimnum(ll))-round(epoch_prestim*fs)+1:cc_events.sample_start(stimnum(ll))+round((epoch_length-epoch_prestim)*fs));
        n=n+1;
    end
end

%% Make figures for specific epoch
elec = 35;
trial = 24;
figure(1)
plot(tt,squeeze(data_epoch(elec,trial,:)))
xlabel('time(s)')
ylabel('amplitude(uV)')
title(sprintf('Electrode %s, stimulating %s',ch{elec},stimchans{trial}))

%% average epochs with same stimulation site

% get the unique number of stimulated pairs:
stimelec = cc_events.electrical_stimulation_site_num(stimevent);
stimnums = cell2mat(cellfun(@str2num, stimelec,'UniformOutput',false));
stimelecs = [sort(stimnums,2) stimcur(stimcur~=0)];
[cc_stimsets,IA,IC] = unique(stimelecs,'rows');

% number of stimuli in each trial
n = histcounts(IC,'BinWidth',0.99);
if any(diff(n) ~=0)
    disp('Not all stimulations are done the same time')
end

%% Averaging epochs
max_stim = max(n);

cc_epoch_sorted = NaN(size(data_epoch,1),max_stim,size(cc_stimsets,1),size(data_epoch,3));

for yy = 1: size(cc_stimsets,1)
    events = find(yy==IC);
    nr_events = length(events);
    
    cc_epoch_sorted(:,1:nr_events,yy,:) = data_epoch(:,events,:);
end
cc_epoch_sorted_avg = squeeze(nanmean(cc_epoch_sorted,2));

%% plot avg epoch
trial = 1;
elec = 3;

figure(1),
plot(tt,squeeze(cc_epoch_sorted(elec,:,trial,:)));
hold on
plot(tt,squeeze(cc_epoch_sorted_avg(elec,trial,:)),'linewidth',2);
hold off
xlabel('time(s)')
ylabel('amplitude(uV)')
title(sprintf('Electrode %s, stimulating %s-%s, %d mA',ch{elec},ch{cc_stimsets(trial,1:2)},cc_stimsets(trial,3)))

%% visually score ERs

ERs = rate_visERssEEG(tt,cc_epoch_sorted,cc_epoch_sorted_avg,ch,cc_stimsets);

%% detect ERs

ERs = detectERssEEG(cc_epoch_sorted_avg,cc_stimsets, fs, epoch_prestim);


