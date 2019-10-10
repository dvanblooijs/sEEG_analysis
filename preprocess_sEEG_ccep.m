%% function for preprocessing SPES in sEEG
% in this function, stimulus pairs with less than 5 stimuli will be removed
% you can choose whether you want to differentiate between
% - negative and positive stimulation (for example: C1-C2 vs C2-C1)
% - different pulse currents (for example: 1 vs 2 mA)

% author: D van Blooijs
% date: Aug 2019

% Sept 2019 - N van Klink - added the possibility to select manually stimulus pairs to keep

function dataBase = preprocess_sEEG_ccep(dataBase,cfg)
epoch_length = cfg.epoch_length;
epoch_prestim = cfg.epoch_prestim;

if exist('minstim','var') == 0
    minstim = 5;
else
    minstim = cfg.minstim;
end

for subj = 1:size(dataBase,2)
    
    fs = dataBase(subj).ccep_header.Fs;
    all_stimcur = str2double(dataBase(subj).tb_events.electrical_stimulation_current)*1000;
    all_stimchans = dataBase(subj).tb_events.electrical_stimulation_site(~isnan(all_stimcur));
    all_stimcur = all_stimcur(~isnan(all_stimcur));
    sample_start = dataBase(subj).tb_events.sample_start(strcmp(dataBase(subj).tb_events.sub_type,'SPES')==1);
    
    %% unique stimulation pairs
    
    % convert stimulus pairs numbers to numbers without unnecessary channels
    stimelecnum = NaN(size(all_stimchans,1),2);
    remove_stim = [];
    for stimp = 1:size(all_stimchans,1)
        stimsplit = strsplit(all_stimchans{stimp},'-');
        for chan = 1:size(stimsplit,2)
            try
                stimelecnum(stimp,chan) = find(strcmp(stimsplit{chan},dataBase(subj).ch)==1);
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
    % use [sort(stimelecnum,2) allstimcur] if you do not want to differentiate direction (positive/negative) of stimulation;
    % use [stimelecnum, allstimcur] if you want to differentiate positive and negative stimulation
    if strcmp(cfg.dir,'yes') && strcmp(cfg.amp,'yes')
        stimelecs = [stimelecnum all_stimcur];
    elseif strcmp(cfg.dir,'yes') && strcmp(cfg.amp,'no')
        stimelecs = stimelecnum;
    elseif strcmp(cfg.dir,'no') && strcmp(cfg.amp,'yes')
        stimelecs = [sort(stimelecnum,2) all_stimcur];
    elseif strcmp(cfg.dir,'no') && strcmp(cfg.amp,'no')
        stimelecs = sort(stimelecnum,2);
    end
    
    [cc_stimsets,~,IC] = unique(stimelecs,'rows');
    
    % check whether all stimulus pairs hase been repeated with the same number stimuli
    n = histcounts(IC,'BinMethod','integers');
    
    if any(diff(n) ~= 0)
        disp('Not all stimulations have been repeated with the same number of stimuli')
        
        x = input('Do you want to manually delete stimulations? [y/n] ','s');
        
        if strcmp(x,'y')
            % check stimpairs and manually remove stimpairs
            avg_stims = round(median(n));
            disp(['Assumed number of stims per pair = ',num2str(avg_stims)]);
            X = find(n~=avg_stims);
            Xstimpnrs = cc_stimsets(X,1:2);
            Xstimpname = dataBase(subj).ch(Xstimpnrs);
            XNstims = n(n~=avg_stims)';
            disp('Check the following stimpairs: ');
            T = table(Xstimpname, XNstims);
            
            % create list of inaccurate stimpairs
            L.stim=NaN;L.name={'none','none'};L.curr=NaN;L.stimnum=NaN;
            for i=1:length(X)
                % find stimulus numbers
                stim = find((stimelecs(:,1)==Xstimpnrs(i,1)) & (stimelecs(:,2)==Xstimpnrs(i,2)));
                L.stim = [L.stim; stim];
                % find stimulus pair names
                name = repmat(Xstimpname(i,:),XNstims(i),1);
                L.name = [L.name; name];
                % find stimulus current
                L.curr = [L.curr; all_stimcur(stim)];
                % find number of stimuli
                stimnum = repmat(XNstims(i,:),XNstims(i),1);
                L.stimnum = [L.stimnum; stimnum];
            end
            
            %concatenate for list view
            Lconc = 'none';
            for i = 2:length(L.stim)
                Lconc= [Lconc; string([num2str(L.stim(i)) '_' L.name{i,1:2} '_' num2str(L.curr(i)) 'mA_' num2str(L.stimnum(i))  ' times stimulated in total'])];
            end
            
            %list to manually select stimulations to keep
            [index,tf]=listdlg('ListString',Lconc,'PromptString','Select stimulations to keep :',...
                'ListSize',[450 500]);
            if tf %only if user made a choice
                P = 1:length(L.stim);
                D = setdiff(P,index); %find index to exclude
                remove_stim_extra = L.stim(D);
                
                %remove from lists
                all_stimchans(remove_stim_extra)=[];
                all_stimcur(remove_stim_extra) = [];
                stimelecs(remove_stim_extra,:) = [];
                sample_start(remove_stim_extra) = [];
                
            end
            
            [cc_stimsets,~,IC] = unique(stimelecs,'rows');
            
            n = histcounts(IC,'BinMethod','integers');
            
        else
            fprintf('No manual selection of stimpairs. Stimpairs with less than %d stimuli are removed.\n',minstim)
            stimremove = find(n<minstim); % remove al stimulation pairs that are stimulated less than 5 times
            
            stimelecs(IC==stimremove,:) = [];
            
            [cc_stimsets,~,IC] = unique(stimelecs,'rows');
            n = histcounts(IC,'BinMethod','integers');
            if any(diff(n) ~= 0)
                fprintf('WARNING: %s some stimulation pairs are stimulated less/more than all others\n',dataBase(subj).sub_label)
            end
            
        end
    end
    
    cc_stimchans = cell(size(cc_stimsets,1),2);
    
    for stimp = 1:size(cc_stimsets,1)
        for chan =1:2
            cc_stimchans{stimp,chan} = dataBase(subj).ch{cc_stimsets(stimp,chan)};
        end
        
    end
    
    max_stim = median(n);
    
    dataBase(subj).cc_stimsets = cc_stimsets;
    dataBase(subj).cc_stimchans = cc_stimchans;
    dataBase(subj).max_stim = max_stim;
    
    stimdif = find(n ~= max_stim);
    if ~isempty(stimdif)
        fprintf('WARNING: %s some stimulation pairs are stimulated less/more than all others\n',dataBase(subj).sub_label)
        
        for stimp =1:size(stimdif,2)
            disp([cc_stimchans{stimdif(stimp),1} '-' cc_stimchans{stimdif(stimp),2}])
        end
    end
    
    %% select epochs
    t = round(epoch_length*dataBase(subj).ccep_header.Fs);
    tt = (1/fs:1/fs:cfg.epoch_length)-cfg.epoch_prestim;

    % allocation
    cc_epoch_sorted = NaN(size(dataBase(subj).data,1),dataBase(subj).max_stim,size(dataBase(subj).cc_stimsets,1),t);
    tt_epoch_sorted = NaN(dataBase(subj).max_stim,size(dataBase(subj).cc_stimsets,1),t); % samplenumbers for each epoch
    
    for elec = 1:size(dataBase(subj).data,1) % for all channels
        for ll = 1:length(dataBase(subj).cc_stimsets) % for all epochs with >4 stimuli
            if strcmp(cfg.dir,'no')
                eventnum1 = find(strcmp(dataBase(subj).tb_events.electrical_stimulation_site,[dataBase(subj).cc_stimchans{ll,1}, '-',dataBase(subj).cc_stimchans{ll,2}]));
                eventnum2 = find(strcmp(dataBase(subj).tb_events.electrical_stimulation_site,[dataBase(subj).cc_stimchans{ll,2}, '-',dataBase(subj).cc_stimchans{ll,1}]));
                eventnum = [eventnum1;eventnum2];
            elseif strcmp(cfg.dir,'yes')
                eventnum = find(strcmp(dataBase(subj).tb_events.electrical_stimulation_site,[dataBase(subj).cc_stimchans{ll,1}, '-',dataBase(subj).cc_stimchans{ll,2}]));
            end
            
            if size(eventnum,1) >dataBase(subj).max_stim
                events = dataBase(subj).max_stim;
            else
                events = size(eventnum,1);
            end
            
            for n=1:events
                cc_epoch_sorted(elec,n,ll,:) = dataBase(subj).data(elec,dataBase(subj).tb_events.sample_start(eventnum(n))-round(epoch_prestim*dataBase(subj).ccep_header.Fs)+1:dataBase(subj).tb_events.sample_start(eventnum(n))+round((epoch_length-epoch_prestim)*dataBase(subj).ccep_header.Fs));
                tt_epoch_sorted(n,ll,:) = dataBase(subj).tb_events.sample_start(eventnum(n))-round(epoch_prestim*dataBase(subj).ccep_header.Fs)+1:dataBase(subj).tb_events.sample_start(eventnum(n))+round((epoch_length-epoch_prestim)*dataBase(subj).ccep_header.Fs);
            end
        end
    end
    
    cc_epoch_sorted_avg = squeeze(nanmean(cc_epoch_sorted,2));
    
    dataBase(subj).cc_epoch_sorted = cc_epoch_sorted;
    dataBase(subj).tt_epoch_sorted = tt_epoch_sorted;
    dataBase(subj).cc_epoch_sorted_avg = cc_epoch_sorted_avg;
    dataBase(subj).tt = tt;
    
    fprintf('...%s has been epoched and averaged... \n',dataBase(subj).sub_label)
    
end