% ER analysis in seeg
% this code is used to load seeg files in BIDS and analyze the early
% responses
% author: D van Blooijs
% date: May 2019

% Aug 2019 - N van Klink - improved script for selecting specific run
% Oct 2019 - D van Blooijs - cleaned code

%% load settings

settings_seeganalysis

%% load specific patient

dataBase = load_sEEGdata(cfg);

%% find unique stimulation pairs and stimulus current and select epochs

dataBase = preprocess_sEEG_ccep(dataBase,cfg);

%% Manually select stimpairs to consider for detection 

dataBase = select_stimpairs(dataBase,cfg);

%% plot avg epoch
subj = 1;
trial = 24;
elec = 3;

epoch_sorted = dataBase(subj).cc_epoch_sorted;
epoch_sorted_avg = dataBase(subj).cc_epoch_sorted_avg;
ch = dataBase(subj).ch;
stimnames = dataBase(subj).cc_stimchans;
if strcmp(cfg.amp, 'yes')
    stimsets = dataBase(subj).cc_stimsets;
else
    stimsets = [dataBase(subj).cc_stimsets NaN(size(dataBase(subj).cc_stimsets,1),1)];
end
tt = dataBase(subj).tt;

figure(1),
plot(tt,squeeze(epoch_sorted(elec,:,trial,:)));
hold on
plot(tt,squeeze(epoch_sorted_avg(elec,trial,:)),'k','linewidth',2);
hold off
xlabel('time(s)')
ylabel('amplitude(uV)')
title(sprintf('Electrode %s, stimulating %s-%s, %1.1f mA',ch{elec},stimnames{trial,:},stimsets(trial,3)))
ylim([-800 800])
xlim([tt(1) tt(end)])

%% visually score ERs
subj=1;
cfg.select = 'yes'; % if a selection of stimulus pairs is made for further analysis

start_rating = 1;
stop_rating = size(dataBase(1).cc_epoch_sorted_avg,2);

for trial=start_rating:stop_rating
    visERs = rate_visERssEEG(dataBase(subj),cfg,trial);
    dataBase(subj).visERs(trial) = visERs;
end


%% detect ERs
subj=1;

dataBase(subj).detERs = detectERssEEG(dataBase(subj),cfg);


