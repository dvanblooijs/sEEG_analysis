% author: Dorien van Blooijs
% May 2019

% this function detects early responses after Single Pulse Electrical Stimulation
% this function should be adapted to be used in sEEG

function ERs = detectERssEEG(pat)

epoch_sorted_avg = pat(1).epoch_sorted_avg;
cc_stimsets = pat(1).cc_stimsets;
fs = pat(1).fs;
epoch_prestim = pat(1).epoch_prestim;

% pre-allocation: variables determined in my thesis
thresh = 2.5;
minSD = 50;
sel = 20;

SDfactor=[];
SD = [];
% extrasamps = 20; % after the stimulus artefact, no data can be read in 20 samples (with 2048Hz)
extrasamps = round(20/2048*fs); %DvB: I'm not sure whether this is correct... should check!

for trial = 1:size(epoch_sorted_avg,2)
    
    detERs = [];
    
    for elec=1:size(epoch_sorted_avg,1) % for each electrode
        
        seep = squeeze(epoch_sorted_avg(elec,trial,:))';
        
        % determine median en sd
        smepmediantotal = median(seep);
        signal_new = seep - smepmediantotal;
        smeprmsbefore = std(signal_new(1,1:round(fs*0.9*epoch_prestim)));
        
        % when sd < minimal sd, than sd = minimal sd
        if smeprmsbefore < minSD
            smeprmsbefore = minSD;
        end
        
        % when the electrode is stimulated
        if ismember(elec,cc_stimsets(trial,1:2))
            sample = 0;
            ampl = 0;
            
            % in other electrode
        elseif ~ismember(elec,cc_stimsets(trial,1:2))
            % find positive and negative peak
            [samppos, amplpos]  = peakfinder(signal_new(1,fs*epoch_prestim+extrasamps:round(fs*(epoch_prestim+0.1))),sel,[],1);
            [sampneg, amplneg] = peakfinder(signal_new(1,fs*epoch_prestim+extrasamps:round(fs*(epoch_prestim + 0.1))),sel,[],-1);
            
            % excluding the first and last sample
            amplpos(samppos==1) = [];
            samppos(samppos==1) =[];
            amplneg(sampneg==1) = [];
            sampneg(sampneg==1) = [];
            amplpos(samppos >= round(fs*0.1)-extrasamps) = [];
            samppos(samppos >= round(fs*0.1)-extrasamps) = [];
            amplneg(sampneg >= round(fs*0.1)-extrasamps) = [];
            sampneg(sampneg >= round(fs*0.1)-extrasamps) = [];
            
            sampleneg =[];
            amplineg =[];
            if ~isempty(sampneg)
                maxamplneg = find(abs(amplneg) == max(abs(amplneg)));
                sampleneg = sampneg(maxamplneg(1))+fs*2+extrasamps;
                amplineg = amplneg(maxamplneg(1));
            end
            
            samplepos =[];
            amplipos =[];
            if ~isempty(samppos)
                maxamplpos = find(abs(amplpos) == max(abs(amplpos)));
                samplepos = samppos(maxamplpos(1)) +fs*2+extrasamps; % the highest and only sample positive
                amplipos = amplpos(maxamplpos(1));           % the highest and only amplitude positive
            end
            
            % determine the amplitude of the largest peak (positive or negative)
            if ~isempty(samplepos) && isempty(sampleneg) % only samppos
                sample = samplepos;
                ampl = amplipos;
            elseif isempty(samplepos) && ~isempty(sampleneg) % only sampneg
                sample = sampleneg;
                ampl = amplineg;
            elseif ~isempty(samplepos) && ~isempty(sampleneg) % both
                if abs(amplipos) > abs(amplineg) % pos>neg
                    sample = samplepos;
                    ampl = amplipos;
                elseif abs(amplineg) > abs(amplipos) % neg>pos
                    sample = sampleneg;
                    ampl = amplineg;
                elseif abs(amplineg) == abs(amplipos) % pos==neg
                    if sampleneg < samplepos % pick location of peak
                        sample = sampleneg;
                        ampl = amplineg;
                    elseif sampleneg > samplepos
                        sample = samplepos;
                        ampl = amplipos;
                    end
                end
            elseif isempty(samplepos) && isempty(sampleneg)
                sample = 0;
                ampl = 0;
            end
        end
        
        % when peak amplitude is saturated, it is deleted
        if abs(ampl) > 3000
            ampl = 0;
        end
        
        % is a peak an ER or not?
        if abs(ampl) > thresh* abs(smeprmsbefore)
            detERs = [detERs elec];
        end
        
        SDfactor(elec) = abs(ampl)/abs(smeprmsbefore);
        SD(elec) = smeprmsbefore;
        Amplitude(elec) = ampl;
    end
    
    ERs(trial).stimpair = cc_stimsets(trial,1:2);
    ERs(trial).stimcur = cc_stimsets(trial,3);
    ERs(trial).detERs = detERs; % ERs for all stimulation pairs
    
end


end


