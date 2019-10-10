% function to read data into dataBase
% author: Dorien van Blooijs
% date: June 2019

function dataBase = load_sEEGdata(cfg)

dataPath = cfg.dataPath;
sub_labels = cfg.sub_labels;
ses_label = cfg.ses_label;
task_label = cfg.task_label;

if isfield(cfg,'run_label')
   run_label = cfg.run_label;
else
    run_label(1:size(sub_labels,2)) = {'run-*'};
end

dataBase = struct([]);
for i=1:size(sub_labels,2)
    D = dir(fullfile(dataPath,sub_labels{i},ses_label,'ieeg',...
        [sub_labels{i} '_' ses_label '_' task_label ,'_',run_label{i}, '_ieeg.eeg']));
    
    % determine run_label
    if ~isfield(cfg,'run_label')
        if size(D,1) == 1
            run_label{i} = D(1).name(strfind(D(1).name,'run-'):strfind(D(1).name,'_ieeg')-1);
        else
            run_label{i} = D(1).name(strfind(D(1).name,'run-'):strfind(D(1).name,'_ieeg')-1);
            fprintf('WARNING: More runs were available for %s_%s_%s, so determine run_label! \n',sub_labels{i},ses_label,task_label)
        end
        
        dataName = fullfile(D(1).folder, D(1).name);
    else
        dataName = fullfile(D(1).folder, D(1).name);
    end
    
    ccep_data = ft_read_data(dataName,'dataformat','brainvision_eeg');
    ccep_header = ft_read_header(dataName);
    
    % load events
    D = dir(fullfile(dataPath,[sub_labels{i}],ses_label,'ieeg',...
        [sub_labels{i} '_' ses_label '_' task_label ,'_',run_label{i},'_events.tsv']));
    
    eventsName = fullfile(D(1).folder, D(1).name);
    
    tb_events = readtable(eventsName,'FileType','text','Delimiter','\t');
    
    % load electrodes
    D = dir(fullfile(dataPath,[sub_labels{i}],ses_label,'ieeg',...
        [sub_labels{i} '_' ses_label ,'_electrodes.tsv']));
    
    elecsName = fullfile(D(1).folder, D(1).name);
    
    tb_electrodes = readtable(elecsName,'FileType','text','Delimiter','\t');
    log_elec_incl = ~strcmp(tb_electrodes.group,'other');
    tb_electrodes = tb_electrodes(log_elec_incl,:);
    
    % load channels
    D = dir(fullfile(dataPath,[sub_labels{i}], ses_label,'ieeg',...
        [sub_labels{i} '_' ses_label '_' task_label ,'_',run_label{i},'_channels.tsv']));
    
    channelsName = fullfile(D(1).folder, D(1).name);
    
    tb_channels = readtable(channelsName,'FileType','text','Delimiter','\t');
    log_ch_incl = strcmp(tb_channels.status_description,'included');
    
    ch = tb_channels.name;
    ch_incl = {ch{log_ch_incl}}';
    
    data = -1*ccep_data(log_ch_incl,:);
    
    dataBase(i).sub_label = sub_labels{i};
    dataBase(i).ses_label = ses_label;
    dataBase(i).task_label = task_label;
    dataBase(i).run_label = run_label{i};
    dataBase(i).dataName = dataName;
    dataBase(i).ccep_header = ccep_header;
    dataBase(i).tb_events = tb_events;
    dataBase(i).tb_channels = tb_channels;
    dataBase(i).tb_electrodes = tb_electrodes;
    dataBase(i).ch = ch_incl;
    dataBase(i).data = data;
    fprintf('...Subject %s has been run...\n',sub_labels{i})
end

disp('All subjects are loaded')
