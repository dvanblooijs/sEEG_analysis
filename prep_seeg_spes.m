% ER analysis in seeg
% this code is used to load seeg files in BIDS and analyze the early
% responses
% author: D van Blooijs, N van Klink
% date: Aug 2019

addpath('/home/nvanklink/Desktop/Github/sEEG_analysis/sEEG_analysis')
addpath('/home/nvanklink/Matlab/fieldtrip-20170212/')
ft_defaults

%% load seeg with SPES from 1 patients

dataPath = '/Fridge/CCEP';
sub_labels = {'RESP0788'};
ses_label = '1';
task_label = 'SPESclin';

i=1;
D = dir(fullfile(dataPath,['sub-' sub_labels{i}],['ses-' ses_label],'ieeg',...
    ['sub-' sub_labels{i} '_ses-' ses_label '_task-' task_label ,'_run-*_ieeg.eeg']));

%if multiple files are available
if size(D,1)>1
    [index,tf] = listdlg('ListString',{D(:).name},'PromptString','Select file to analyze:',...
        'ListSize',[500 300],'SelectionMode','single');
    dataFName = fullfile(D(index).folder, D(index).name);
else
    dataFName = fullfile(D(1).folder, D(1).name);
end

ccep_dataraw = ft_read_data(dataFName,'dataformat','brainvision_eeg');

ccep_header = ft_read_header(dataFName);
fs = ccep_header.Fs;

% load events
D = dir(fullfile(dataPath,['sub-' sub_labels{i}],['ses-' ses_label],'ieeg',...
    ['sub-' sub_labels{i} '_ses-' ses_label '_task-' task_label ,'_run-*_events.tsv']));

if size(D,1)>1
    eventsFName = fullfile(D(index).folder, D(index).name);
else
    eventsFName = fullfile(D(1).folder, D(1).name);
end

tb_events = readtable(eventsFName,'FileType','text','Delimiter','\t');

% load channels
D = dir(fullfile(dataPath,['sub-' sub_labels{i}],['ses-' ses_label],'ieeg',...
    ['sub-' sub_labels{i} '_ses-' ses_label '_task-' task_label ,'_run-*_channels.tsv']));

if size(D,1)>1
    channelsFName = fullfile(D(index).folder, D(index).name);
else
    channelsFName = fullfile(D(1).folder, D(1).name);
end

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
% including option to manually exclude stimulations with different number
% of repetitions

all_stimcur = str2double(tb_events.electrical_stimulation_current)*1000;
all_stimchans = tb_events.electrical_stimulation_site(~isnan(all_stimcur));
all_stimcur = all_stimcur(~isnan(all_stimcur));
sample_start = pat(1).sample_start;

% convert stimulus pairs numbers to numbers without unnecessary channels
stimelecnum = NaN(size(all_stimchans,1),2);
remove_stim = [];
for stimp =1:size(all_stimchans,1)
    stimsplit = strsplit(all_stimchans{stimp},'-');
    for elec = 1:size(stimsplit,2)
        try
            stimelecnum(stimp,elec) = find(strcmp(pat(1).ch,stimsplit{elec}));
        catch % if stimp elec is not in channel list (e.g. is bad channel)
            remove_stim = [remove_stim, stimp];
        end
    end
end

%remove stimpairs with bad channels
all_stimchans(remove_stim)=[];
all_stimcur(remove_stim) = [];
stimelecnum(remove_stim,:) = [];
sample_start(remove_stim) = [];

% get the unique number of stimulated pairs:
 % use [sort(stimelecnum,2) all_stimcur] if you do not want to differentiate direction (positive/negative) of stimulation; 
 % use [stimelecnum, all_stimcur] if you want to differentiate positive and negative stimulation
p=input('Pos and Neg evaluation separate? [y/n] ','s');
if strcmp(p,'y') 
    stimelecs = [stimelecnum, all_stimcur];
else
    stimelecs = [sort(stimelecnum,2) all_stimcur];
    for i=1:size(stimelecs,1)
        all_stimchans{i} = strcat(pat.ch(stimelecs(i,1)), '-', pat.ch(stimelecs(i,2)));
    end
end
[cc_stimsets,IA,IC] = unique(stimelecs,'rows');

% number of stimuli in each trial
n = histcounts(IC,'BinMethod','integers');

if any(diff(n) ~=0)
    disp('Not all stimulations have been repeated with the same number of stimuli')
    
    % check stimpairs and manually remove stimpairs
    avg_stims = round(mean(n));
    disp(['Assumed number of stims per pair = ',num2str(avg_stims)]);
    X = find(n~=avg_stims);
    Xstimpnrs = cc_stimsets(X,1:2);
    Xstimpname = pat(1).ch(Xstimpnrs);
    XNstims = n(n~=avg_stims)';
    disp('Check the following stimpairs: ');
    T = table(Xstimpname, XNstims)
        
    x = input('Do you want to manually delete stimulations? [y/n] ','s');
    
    if strcmp(x,'y')
        %create list of inaccurate stimpairs
        L.stim=[];L.name=[];L.curr=[];
        for i=1:length(X)
            stim = find((stimelecs(:,1)==Xstimpnrs(i,1)) & (stimelecs(:,2)==Xstimpnrs(i,2)));
            L.stim = [L.stim; stim];
            name = repmat([Xstimpname(i,:)],XNstims(i),1);
            L.name = [L.name; name];
            L.curr = [L.curr; all_stimcur(stim)];
        end

        %concatenate for list view
        Lconc = [];
        for i = 1:length(L.stim)
            Lconc= [Lconc; string([num2str(L.stim(i)) '_' L.name{i,1:2} '_' num2str(L.curr(i)) 'mA'])];
        end

        %list to manually select stimulations to keep
        [index,tf]=listdlg('ListString',Lconc,'PromptString','Select stimulations to keep:',...
            'ListSize',[200 500]);
        if tf %only if user made a choice
            P = 1:length(L.stim);
            D = setdiff(P,index); %find index to exclude
            remove_stim_extra = L.stim(D);

            %remove from lists
            all_stimchans(remove_stim_extra)=[];
            all_stimcur(remove_stim_extra) = [];
            stimelecnum(remove_stim_extra,:) = [];
            sample_start(remove_stim_extra) = [];

            %recompute summary list  
            % use [sort(stimelecnum,2) all_stimcur] if you do not want to differentiate direction (positive/negative) of stimulation; 
            % use [stimelecnum, all_stimcur] if you want to differentiate positive and negative stimulation
            if strcmp(p,'y') 
                stimelecs = [stimelecnum, all_stimcur];
            else 
                stimelecs = [sort(stimelecnum,2) all_stimcur];
            end
            [cc_stimsets,IA,IC] = unique(stimelecs,'rows');

            n = histcounts(IC,'BinMethod','integers');
        end
    end
end

pat(1).all_stimchans = all_stimchans;
pat(1).all_stimcur = all_stimcur;
pat(1).cc_stimsets = cc_stimsets;
pat(1).stimnum_max = max(n);
pat(1).sample_start = sample_start;
pat(1).all_stimnrs = stimelecs(:,1:2);
pat(1).IC = IC;

% remove all irrelevant variables
clearvars -except pat tb_channels

%% Manually select stimpairs to consider for detection
% create list of stimpairs
conc=[];
stimpname = pat(1).ch(pat(1).cc_stimsets(:,1:2));
for i=1:size(pat(1).cc_stimsets,1)
conc = [conc; string([stimpname{i,1:2} '_' num2str(pat(1).cc_stimsets(i,3)) 'mA'])];
end

[index,tf]=listdlg('ListString',conc,'PromptString','Select stimulations to keep:',...
    'ListSize',[200 500]);
if tf %only if user made a choice
    keep_stimsets = pat(1).cc_stimsets(index,:); %find stimpairs
    
    keep_stims = [];
    for i=1:length(index) %find individual stims
        keep_stims = [keep_stims; find((pat(1).all_stimnrs(:,1)==keep_stimsets(i,1)) & (pat(1).all_stimnrs(:,2)==keep_stimsets(i,2)))];
    end
    
    %reduce patient variable
    pat(1).all_stimchans = pat(1).all_stimchans(keep_stims);
    pat(1).all_stimcur = pat(1).all_stimcur(keep_stims);
    pat(1).cc_stimsets = keep_stimsets;
    pat(1).sample_start = pat(1).sample_start(keep_stims);
    pat(1).all_stimnrs = pat(1).all_stimnrs(keep_stims,:);
    [~,~,IC] = unique(pat(1).all_stimnrs,'rows');
    pat(1).IC = IC;
end

% remove all irrelevant variables
clearvars -except pat tb_channels

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

clearvars -except pat tb_channels


%% Make figures for specific epoch
elec = 35;
trial = 3; 
data_epoch = pat(1).data_epoch;
ch_incl = pat(1).ch;
all_stimchans = pat(1).all_stimchans;
tt= pat(1).tt;

figure(1)
plot(tt,squeeze(data_epoch(elec,trial,:)))
xlabel('time(s)')
ylabel('amplitude(uV)')
title(sprintf('Electrode %s, stimulating %s, %d mA',ch_incl{elec},string(all_stimchans{trial}),pat(1).all_stimcur(trial)))
ylim([-1600 1600])

clearvars -except pat tb_channels


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

clearvars -except pat tb_channels


%% plot avg epoch
trial = 1;
elec = 11;

epoch_sorted = pat(1).epoch_sorted;
epoch_sorted_avg = pat(1).epoch_sorted_avg;
ch = pat(1).ch;
stimnames = pat(1).cc_stimchans;
stimsets = pat(1).cc_stimsets;
tt = pat(1).tt;

figure,
plot(tt,squeeze(epoch_sorted(elec,:,trial,:)));
hold on
plot(tt,squeeze(epoch_sorted_avg(elec,trial,:)),'k','linewidth',2);
hold off
xlabel('time(s)')
ylabel('amplitude(uV)')
title(sprintf('Electrode %s, stimulating %s, %1.1f mA',ch{elec},string(stimnames{trial}),stimsets(trial,3)))
ylim([-800 800])

clearvars -except pat tb_channels

%% visually score ERs

start_rating = 1;
stop_rating = size(pat(1).epoch_sorted_avg,2);

for trial=start_rating:stop_rating
    visERs = rate_visERssEEG(pat,trial);
    pat(1).visERs(trial) = visERs;
end

visERs = pat.visERs;
% clearvars -except pat tb_channels


%% detect ERs

pat(1).detERs = detectERssEEG(pat);

clearvars -except pat tb_channels

