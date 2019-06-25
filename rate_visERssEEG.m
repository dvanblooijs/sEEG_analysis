% author: Dorien van Blooijs
% date: May 2019

% this function scores ERs. You should press 'y' when an ER is observed.
% When there is no ER, you can press enter.
% with start_rating and stop_rating, you can choose where to start with
% rating. It is recommended to choose parts of the data, so that it saves
% after you finishes each part. 


function ERs = rate_visERssEEG(pat,start_rating,stop_rating)

tt = pat(1).tt;
epoch_sorted = pat(1).epoch_sorted;
epoch_sorted_avg = pat(1).epoch_sorted_avg;
ch = pat(1).ch;
cc_stimsets = pat(1).cc_stimsets;
cc_stimchans = pat(1).cc_stimchans;
stimcur = pat(1).cc_stimsets(:,3);

ERs = struct;
% for each stimulus
for trial=start_rating:stop_rating
    ER = [];
    for elec =1:size(epoch_sorted_avg,1)
        if ~ismember(elec,cc_stimsets(trial,1:2))
            % figure with left the epoch, and right zoomed in
            H=figure(1);
            H.Units = 'normalized';
            H.Position = [0.13 0.11 0.77 0.8];
            
            subplot(1,2,1)
            plot(tt,squeeze(epoch_sorted(elec,:,trial,:)),'r');
            hold on
            plot(tt,squeeze(epoch_sorted_avg(elec,trial,:)),'k','linewidth',2);
            hold off
            xlim([-2 2])
            ylim([-2000 2000])
            xlabel('time(s)')
            ylabel('amplitude(uV)')
            title(sprintf('Electrode %s, stimulating %s, %1.1f mA',ch{elec},cc_stimchans{trial},stimcur(trial)))
            
            subplot(1,2,2)
            plot(tt,squeeze(epoch_sorted(elec,:,trial,:)),'r');
            hold on
            plot(tt,squeeze(epoch_sorted_avg(elec,trial,:)),'k','linewidth',2);
            hold off
            xlim([-0.5 1])
            ylim([-750 750])
            title('Zoomed average signal')
            xlabel('Time (s)')
            ylabel('Voltage (uV)')
            x = input('ER? [y/n] ','s');
            if strcmp(x,'y')
                ER = [ER elec];
            end
        end
    end
    ERs(trial).stimpair = cc_stimsets(trial,1:2);
    ERs(trial).vis = ER;
    ERs(trial).stimcur = cc_stimsets(trial,3);
    
end
end


