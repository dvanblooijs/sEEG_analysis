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
sub_labels = {'RESP0773'};
ses_label = '1';
task_label = 'SPESclin';

i=1;
D = dir(fullfile(dataPath,['sub-' sub_labels{i}],['ses-' ses_label],'ieeg',...
    ['sub-' sub_labels{i} '_ses-' ses_label '_task-' task_label ,'_run-*_ieeg.eeg']));

dataFName = fullfile(D(1).folder, D(1).name);

ccep_dataraw = ft_read_data(dataFName,'dataformat','brainvision_eeg');

ccep_header = ft_read_header(dataFName);
fs = ccep_header.Fs;

% load events
D = dir(fullfile(dataPath,['sub-' sub_labels{i}],['ses-' ses_label],'ieeg',...
    ['sub-' sub_labels{i} '_ses-' ses_label '_task-' task_label ,'_run-*_events.tsv']));

eventsFName = fullfile(D(1).folder, D(1).name);

tb_events = readtable(eventsFName,'FileType','text','Delimiter','\t');

% load channels
D = dir(fullfile(dataPath,['sub-' sub_labels{i}],['ses-' ses_label],'ieeg',...
    ['sub-' sub_labels{i} '_ses-' ses_label '_task-' task_label ,'_run-*_channels.tsv']));

channelsFName = fullfile(D(1).folder, D(1).name);

tb_channels = readtable(channelsFName,'FileType','text','Delimiter','\t');
% ch = tb_channels.name;

% remove all unnecessary electrodes
log_ch_incl = strcmp(tb_channels.status_description,'included');
ch_incl = tb_channels.name(log_ch_incl);

ccep_data = -1*ccep_dataraw(log_ch_incl,:); %*-1 because that's how the ECoG is displayed and visualized
% put all information to use in one struct called "pat"
pat(i).dataFName = dataFName;
pat(i).RESPnum = sub_labels{:};
pat(i).ch = ch_incl;
pat(i).fs = fs;
pat(i).ccep_data = ccep_data;
pat(i).sample_start = tb_events.sample_start(strcmp(tb_events.sub_type,'SPES')==1);

% remove all irrelevant variables
clearvars -except pat tb_events tb_channels

%% find unique stimulation pairs and stimulus current

all_stimcur = str2double(tb_events.electrical_stimulation_current)*1000;
all_stimchans = tb_events.electrical_stimulation_site(~isnan(all_stimcur));
all_stimcur = all_stimcur(~isnan(all_stimcur));

% convert stimulus pairs numbers to numbers without unnecessary channels
stimelecnum = NaN(size(all_stimchans,1),2);
for stimp =1:size(all_stimchans,1)
    stimsplit = strsplit(all_stimchans{stimp},'-');
    for elec = 1:size(stimsplit,2)
       stimelecnum(stimp,elec) = find(strcmp(pat(1).ch,stimsplit{elec}));
    end
end

% get the unique number of stimulated pairs:
 % use [sort(stimelecnum,2) allstimcur] if you do not want to differentiate direction (positive/negative) of stimulation; 
 % use [stimelecnum, allstimcur] if you want to differentiate positive and negative stimulation
stimelecs = [stimelecnum, all_stimcur];
[cc_stimsets,IA,IC] = unique(stimelecs,'rows');

% number of stimuli in each trial
n = histcounts(IC,'BinMethod','integers');
if any(diff(n) ~=0)
    disp('Not all stimulations have been repeated with the same number of stimuli')
end

pat(1).all_stimchans = all_stimchans;
pat(1).all_stimcur = all_stimcur;
pat(1).cc_stimsets = cc_stimsets;
pat(1).IC = IC;
pat(1).stimnum_max = max(n);

% remove all irrelevant variables
clearvars -except pat tb_events tb_channels

%% epoch files in 2spre-2spost stimulus 
epoch_length = 4; % in seconds, -2:2
epoch_prestim = 2;
fs = pat(1).fs;
data_epoch = zeros(size(pat(1).ch,1),size(pat(1).all_stimchans,1),round(epoch_length*fs)); % preallocation: [number of stimuli x epoch length]

tt = (1/fs:1/fs:epoch_length) - epoch_prestim;

for elec = 1:size(pat(1).ch,1) % for all channels
    for ll = 1:size(pat(1).all_stimchans,1) % for all single stimulations
        data_epoch(elec,ll,:) = pat(1).ccep_data(elec,pat(1).sample_start(ll)-round(epoch_prestim*fs)+1:pat(1).sample_start(ll)+round((epoch_length-epoch_prestim)*fs));                
    end
end

pat(1).epoch_prestim = epoch_prestim;
pat(1).data_epoch = data_epoch;
pat(1).tt = tt;

clearvars -except pat tb_events tb_channels


%% Make figures for specific epoch
elec = 35;
trial = 24; % this is in systemplus also the triggernumber-1000
data_epoch = pat(1).data_epoch;
ch_incl = pat(1).ch;
all_stimchans = pat(1).all_stimchans;
tt= pat(1).tt;

figure(1)
plot(tt,squeeze(data_epoch(elec,trial,:)))
xlabel('time(s)')
ylabel('amplitude(uV)')
title(sprintf('Electrode %s, stimulating %s, %d mA',ch_incl{elec},all_stimchans{trial},pat(1).all_stimcur(trial)))
ylim([-1600 1600])

clearvars -except pat tb_events tb_channels


%% Averaging epochs

cc_epoch_sorted = NaN(size(pat(1).data_epoch,1),pat(1).stimnum_max,size(pat(1).cc_stimsets,1),size(pat(1).data_epoch,3));
cc_stimchans = cell(size(pat(1).cc_stimsets,1),1);

data_epoch = pat(1).data_epoch;

for yy = 1: size(pat(1).cc_stimsets,1)
    events = find(yy==pat(1).IC); % find where the sorted stim channels are really stimulated
    nr_events = length(events);
    
    cc_epoch_sorted(:,1:nr_events,yy,:) = data_epoch(:,events,:); % [channels,stimulus1-10,trials,time]
    cc_stimchans{yy,1} = pat(1).all_stimchans{events(1)};
end
cc_epoch_sorted_avg = squeeze(nanmean(cc_epoch_sorted,2));

pat(1).epoch_sorted = cc_epoch_sorted;
pat(1).epoch_sorted_avg = cc_epoch_sorted_avg;
pat(1).cc_stimchans = cc_stimchans;

clearvars -except pat tb_events tb_channels


%% plot avg epoch
trial = 54;
elec = 3;

epoch_sorted = pat(1).epoch_sorted;
epoch_sorted_avg = pat(1).epoch_sorted_avg;
ch = pat(1).ch;
stimnames = pat(1).cc_stimchans;
stimsets = pat(1).cc_stimsets;
tt = pat(1).tt;

figure(1),
plot(tt,squeeze(epoch_sorted(elec,:,trial,:)));
hold on
plot(tt,squeeze(epoch_sorted_avg(elec,trial,:)),'k','linewidth',2);
hold off
xlabel('time(s)')
ylabel('amplitude(uV)')
title(sprintf('Electrode %s, stimulating %s, %1.1f mA',ch{elec},stimnames{trial},stimsets(trial,3)))
ylim([-800 800])

clearvars -except pat tb_events tb_channels

%% visually score ERs

start_rating = 77;%1;
stop_rating = size(pat(1).epoch_sorted_avg,2);

for trial=start_rating:stop_rating
    visERs = rate_visERssEEG(pat,trial);
    pat(1).visERs(trial) = visERs;
end

% clearvars -except pat tb_events tb_channels


%% detect ERs

pat(1).detERs = detectERssEEG(pat);

clearvars -except pat tb_events tb_channels

