% author: Dorien van Blooijs
% May 2019

% this function scores ERs. You should press 'y' when an ER is observed.


function ERs = rate_visERssEEG(tt,cc_epoch_sorted,cc_epoch_sorted_avg,ch,cc_stimsets)

ERs = struct;
% for each stimulus
for trial=1:size(cc_epoch_sorted_avg,2)
    ER = [];
    for elec =1:size(cc_epoch_sorted_avg,1)
        if ~ismember(elec,cc_stimsets(trial,1:2))
            % figure with left the epoch, and right zoomed in
            H=figure(1);
            H.Units = 'normalized';
            H.Position = [0.13 0.11 0.77 0.8];
            
            subplot(1,2,1)
            plot(tt,squeeze(cc_epoch_sorted(elec,:,trial,:)),'r');
            hold on
            plot(tt,squeeze(cc_epoch_sorted_avg(elec,trial,:)),'k','linewidth',2);
            hold off
            xlim([-2 2])
            ylim([-2000 2000])
            xlabel('time(s)')
            ylabel('amplitude(uV)')
            title(sprintf('Electrode %s, stimulating %s-%s, %d mA',ch{elec},ch{cc_stimsets(trial,1:2)},cc_stimsets(trial,3)))
            
            subplot(1,2,2)
            plot(tt,squeeze(cc_epoch_sorted(elec,:,trial,:)),'r');
            hold on
            plot(tt,squeeze(cc_epoch_sorted_avg(elec,trial,:)),'k','linewidth',2);
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
end
end


