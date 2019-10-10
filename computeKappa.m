% compute kappa for two observers
% for visual SPES-ER analysis

V1 = visERs_N;
V2 = visERs_D;

clear visERs

%create empty scoring matrix
score1 = zeros(size(pat.epoch_sorted_avg,1),size(pat.epoch_sorted_avg,2));
score2 = score1;

%check if all stimpairs are in the V1 sturcture
if sum(~cellfun(@isempty,{V1.stimpair})) == size(pat.epoch_sorted_avg,2) && ...
        sum(~cellfun(@isempty,{V2.stimpair})) == size(pat.epoch_sorted_avg,2)

    %include scores in scoring matrix
    for i=1:sum(~cellfun(@isempty,{V1.stimpair}))
        %run through stimpairs
        score1(V1(i).vis,i) = 1;
        score2(V2(i).vis,i) = 1;
    end

    % get 
    agreeY = (score1==1 & score2==1);
    agreeN = (score1==0 & score2==0);
    Y1N2 = (score1==1 & score2==0);
    Y2N1 = (score1==0 & score2==1);

    TagreeY = sum(agreeY(:));
    TagreeN = sum(agreeN(:));
    TY1N2 = sum(Y1N2(:));
    TY2N1 = sum(Y2N1(:));

    x(1,1) = TagreeY;
    x(1,2) = TY1N2;
    x(2,1) = TY2N1;
    x(2,2) = TagreeN;

    %compute kappa
    confLevel= 0.95;
    dispResults =0; %Whether to show results in stdoutput
    infoKAPPA = kappa(x,0,1-confLevel,dispResults);

else
    error('number of stimpairs in Vis_ERs is not correct');
end

%% put effort into disgreement 
% check disagreement matrices
mtx = Y1N2 | Y2N1;

for i=1:size(mtx,2)
    trial = i;
    ER = [];
    elec_all = find(mtx(:,i));

    for r = 1:length(elec_all)
        elec = elec_all(r);

    figure,
    plot(pat.tt,squeeze(pat.epoch_sorted(elec,:,trial,:)));
    hold on
    plot(pat.tt,squeeze(pat.epoch_sorted_avg(elec,trial,:)),'k','linewidth',2);
    hold off
    xlabel('time(s)')
    ylabel('amplitude(uV)')
    title(sprintf('Electrode %s, stimulating %s, %1.1f mA',pat.ch{elec},pat.cc_stimchans{trial},pat.cc_stimsets(trial,3)))
    ylim([-800 800])

    x = input('ER? [y/n] ','s');
    if strcmp(x,'y')
        mtx2(elec,trial) = 1;
        ER = [ER elec];
    end
    end
    close all;
    %reconstruct scored ERs after visual check
    %add already agreed ERs    
    visERs(trial).stimpair = pat.cc_stimsets(trial,1:2);
    visERs(trial).vis = find(agreeY(:,trial));
    visERs(trial).stimcur = pat.cc_stimsets(trial,3);    

    %add newly checked ERs
    visERs(trial).vis = sort([visERs(trial).vis; ER']);
end