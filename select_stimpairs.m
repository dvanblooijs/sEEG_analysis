% function to select stimulus pairs for further analysis
% author: N van Klink
% date: Aug 2019

% Oct 2019 - D van Blooijs - cleaned up code and put it in function

function dataBase = select_stimpairs(dataBase, cfg)

for subj = 1:size(dataBase,2)
    
    % create list of stimpairs
    conc=[];
    stimpname = dataBase(subj).ch(dataBase(subj).cc_stimsets(:,1:2));
    
    if strcmp(cfg.amp,'yes')
        for i=1:size(dataBase(subj).cc_stimsets,1)
            conc = [conc; string([stimpname{i,1:2} '_' num2str(dataBase(subj).cc_stimsets(i,3)) 'mA'])];
        end
    else
        for i=1:size(dataBase(subj).cc_stimsets,1)
            conc = [conc; string([stimpname{i,1:2}])];
        end
    end
    
    [index,tf]=listdlg('ListString',conc,'PromptString','Select stimulations to keep:',...
        'ListSize',[200 500]);
    if tf %only if user made a choice
        
        %reduce patient variable
        dataBase(subj).select_stimchans = dataBase(subj).cc_stimchans(index,:);
        dataBase(subj).select_stimsets = dataBase(subj).cc_stimsets(index,:);
        
        % select data
        dataBase(subj).select_epoch_sorted = dataBase(subj).cc_epoch_sorted(index,:,:,:);
        dataBase(subj).select_epoch_sorted_avg = squeeze(mean(dataBase(subj).select_epoch_sorted,2));
        
    else
        disp('ERROR: no stimulus pairs were selected')
    end
end