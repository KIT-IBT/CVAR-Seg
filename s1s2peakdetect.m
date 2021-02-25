% -------------------------------------------------------
%
%    s1s2peakdetect  - CVAR-Seg main class
%
%    Ver. 1.0.0
%
%    Created:           Mark Nothstein (25.02.2020)
%    Last modified:     Mark Nothstein (25.02.2020)
%
%    Institute of Biomedical Engineering
%    Karlsruhe Institute of Technology
%
%    http://www.ibt.kit.edu
%
%    Copyright 2000-2020 - All rights reserved.
%
% ------------------------------------------------------
%
% location(namein,sampleratein,xyzin,elin,cathin,timein)
% class storing geometric data. 
%
% Expected input are unipolar positions of electrodes.
% Functionality to calc estimated bipolar positions, distance between
% catheter electrodes (uni or bipolar)
%
%
% Inputs:
% Required:
%      'signal_in_bi'   array - (samples x value)
%      'samplerate'     int - (samplerate Hz, e.g 1000)
%      'leads_in'       cellarray ({'CS';'CS';...})
%      'elecs_in'       cellarray ({'1-2';'2-3';...})
%      'type'           'bipolar','unipolar' (legacy -> only bipolar is supported)
%      'location'       location class -> see classes folder
%      'map'            map class -> see classes folder
%      'ECG'            location class -> see classes folder
% Optional:
%      'signal_in_uni'  array - (samples x value)
%      'leads_in_uni'   cellarray ({'CS';'CS';...})
%      'elecs_in_uni'   cellarray ({'1';'2';...})
%      'name'           string 
%      's1time'         double
%      's2times'        array double
%      'ERP'            double
%The values s1time,s2time,ERP are never used and are only ment as a way to
%store all measurement data in one class
%
% Outputs:
%           location class
%
%
% Example Usage:
%    see run_cvar_bench
%
% Revision history:

classdef s1s2peakdetect < handle
    properties
        signal       %struct storing signal_in_bi,bi etc..
        
        
        signalflag
        isfiltered
        samplerate
        leads
        elecs
        locations
        map
        ecg
        
        stimchan
        timeline
        
        name
        init
        weight
        threshold
        guess_siglength_stim
        variables
        
        classcheck
        
        S1_time
        numof_stimuli_in_one_block
        numof_S1_in_one_block
        S1_blocknum
        S1_blocknum_locs
        s2_steps
        s2_fs
        erp_id
        bad_aa_blocks
        
        cv_plot
        amp_plot
        
        atrial_segments
        signal_properties
        nonstimchan_signal_properties
        
        %manual method
        manual_segments_stimchan_list % not tested or developed further
        settings
    end
    
    methods
        function obj = s1s2peakdetect(signal_in_bi,samplerate,leads_in,elecs_in,type,varargin)
            %Constructs
            p = inputParser;
            addRequired(p,'signal_in_bi',@(x) ~isempty(x) && isnumeric(x));
            addRequired(p,'samplerate',@(x) ~isempty(x) && isnumeric(x) && size(x,1)==1 && size(x,2)==1);
            addRequired(p,'leads_in',@(x) ~isempty(x) && max(size(x))==min(size(signal_in_bi)) );
            addRequired(p,'elecs_in',@(x) ~isempty(x) && max(size(x))==min(size(signal_in_bi)) );
            addRequired(p,'type',@(x) ~isempty(x) && ischar(x));
            addParameter(p,'signal_in_uni','',@(x) ~isempty(x) && isnumeric(x));
            addParameter(p,'leads_in_uni',@(x) ~isempty(x) && max(size(x))==min(size(signal_in_uni)) );
            addParameter(p,'elecs_in_uni',@(x) ~isempty(x) && max(size(x))==min(size(signal_in_uni)) );
            addParameter(p,'name','',@(x) ischar(x) && ~isempty(x) );
            addParameter(p,'s1time','',@(x) isnumeric(x) && ~isempty(x) );
            addParameter(p,'s2times','',@(x) isnumeric(x) && ~isempty(x) );
            addParameter(p,'ERP','',@(x) isnumeric(x) && ~isempty(x) && size(x,1)==1 && size(x,2)==1 );
            addParameter(p,'location','',@(x) isa(x,'location') && ~isempty(x) );
            addParameter(p,'map','',@(x) isa(x,'map') && ~isempty(x) );
            addParameter(p,'ECG','',@(x) ~isempty(x) );
            try
                p.parse( signal_in_bi,samplerate,leads_in,elecs_in,type,varargin{:} );
            catch MyError
                rethrow(MyError);
            end
            
            %needed values
            obj.signal.in_bi = p.Results.signal_in_bi;
            obj.samplerate = p.Results.samplerate;
            obj.leads.in = p.Results.leads_in;
            obj.elecs.in = p.Results.elecs_in;
            obj.signalflag = p.Results.type;
            obj.locations = p.Results.location;
            
            %optional values
            obj.signal.in_uni = p.Results.signal_in_uni;
            obj.leads.in_uni = p.Results.leads_in_uni;
            obj.elecs.in_uni = p.Results.elecs_in_uni;
            obj.name = p.Results.name;
            obj.classcheck.ERP = p.Results.ERP;
            obj.classcheck.s1time = p.Results.s1time;
            %s2
            if ~isempty(p.Results.s2times) && size(p.Results.s2times,1)<size(p.Results.s2times,2) %s2times given
                obj.classcheck.s2times = p.Results.s2times';
            elseif ~isempty(p.Results.s2times)
                obj.classcheck.s2times = p.Results.s2times;
            else
                obj.classcheck.s2times = p.Results.s2times;
            end
            obj.map = p.Results.map;
            obj.ecg.in = p.Results.ECG;
            clear p varargin signal.in samplerate isunipolar leads_in
            
            %settings
            obj.settings.manualmode = 0;
            obj.settings.ploton = 0;
            if obj.settings.manualmode %since no check possible without plots
                obj.settings.ploton = 1; 
            end
            
            %initialization   
            obj.init.k = 0.01*obj.samplerate/1000;                  %high k means less peaks
            obj.init.step_k = 0.01;
            
            obj.init.lp_uni = 500;                                  %change this once using unipolar
            obj.init.hp_uni = 5;
            
            obj.init.lp_bi = 250;
            obj.init.hp_bi = 1;
            
            obj.init.filt_order = 2;
            
            obj.init.min_CV = 10;                                   %units: mm/s
            obj.init.s = 20;                                        %units: mm
            
            obj.init.highhp = 450;                                  %used to find stimulipeaks, should be larger than 300! & lower than whatever filter set in init.lp
            
            obj.init.S1_guess_good = false;
            obj.init.S1_block_guess_good = false;                   
            obj.init.foundS1_sequence = false;
            obj.init.s1s2found = false;
            obj.init.step_weighting_increase = 1;                   %dependant on your weighting scheme,for mine should be int!
            
            
            %weighting factors
            obj.weight.detected_in_stim_chan = 2;                   %detected peaks in stimchan most likely primarily stimulations
            obj.weight.detected_in_chan = 1;
            obj.weight.hp_detected_in_stim_chan = 20;
            obj.weight.hp_detected_in_chan = 10;
            obj.weight.dist = 100;
            obj.weight.diff = 100;
            obj.weight.amp = 100;
            
            %thresholds
            obj.threshold.deviation_S1_guess = 10;                  %units: samples, by how many samples the neighbour peaks are alowed to differ
            obj.threshold.max_itterations_S1_guess = 100;
            obj.threshold.max_neighbour = 5;                        %number of multiples of distancees that should be evaluated by histogram,e.g. 3*s1dist
            obj.threshold.s1_block_guess_variance = 2;              %programm guesses how many stimuli are in an s1 train. To test this we add and subtract a variance and test if result gets better or worse. eg guess is 6, variance is 3 we test (3 4 5 6 7 8 9)
            obj.threshold.max_deviation_S1_block_guess = 100;       %units: samples
            obj.threshold.deviation_S1_block_guess = 0.08*obj.samplerate;      %by how many samples the S1 guess is allowed to be wrong
            obj.threshold.increase_S1_block_deviation= 2;           %increments (units: samples) that the threshold should be lowered
            obj.threshold.stimuluslength_fs = (0.002+0.015)*obj.samplerate; %most stimulus generators have an exponential capacitor decay of ~6ms (UHS30) or 15ms(UHS3000)
            obj.threshold.stimuluslength_fs_deviation = round(obj.threshold.stimuluslength_fs * 0.5);
            obj.threshold.ERP_ratio_perc = 0.2;
            
            %variables changing
            obj.variables.k = obj.init.k;
            obj.variables.method_s1s2detect = 'wavelet';
            %obj.start_detection
            
        end
        
                
        function start_detection(obj)
            %starts detection with unfiltered bipolar navx signals
            
            %baselineremoval
            %===============
            obj.getbipolarsignalsandleads();
            obj.splitcatheters();
            obj.getstimchannel();
            if obj.settings.manualmode==1
                obj.presegmentsignal();
            else
                obj.ecg.filt.ydata = obj.ecg.in.ydata;
            end 
            obj.detectEKGpeaks();
            
            %remove baselinedrift etc... Extreme filtering is done later on
            %in each function depending on what it focusses on.
            obj.prefilter_signals();
            
            obj.guessS1time();
            if obj.settings.ploton
                obj.timeline.plotsupersegments(obj.signal.bi_stimchan(:,obj.stimchan.distant_el_id),0)
            end
            
            obj.high_lp_filt_weighting(obj.variables.k);
            if obj.settings.ploton
                obj.timeline.plotsupersegmentweights(obj.signal.bi_stimchan);
            end
            
            obj.guessblocklocations
            obj.deletesegmentsoutsideblocks(); %dependant on perfect block detection!
            
            obj.adddistanceweighting();
            obj.addthresholdweighting;
            obj.adddistanceweighting_old(); 
            obj.findS1S2peaks();
            
            if strcmp(obj.variables.method_s1s2detect,'wavelet')
                obj.refine_segments();
            end
            obj.gets2times();
            obj.checkVF()
            obj.getPwavestats();
            obj.decideStimTemplate();
            obj.gets1activity();
            obj.getatrialactivity(); %something goes wrong here!!!
            obj.combine_segmented_sig();
            obj.getatrialamplitudes();
            obj.gets1activity_nonstimchan();
            obj.getatrialactivity_nonstimchan();
            obj.getatrialamplitudes_nonstimchan();
            obj.detectERP(obj.threshold.ERP_ratio_perc);
            
            %unipolar part
            if ~isempty(obj.signal.in_uni)
                obj.getatrialactivity_unipolar();
            end
            
            try
            obj.getlocations();
            catch %Here for now cause of benchmarkpaper
            end
            %obj.plotLassoCV_theo(); %not really needed/necessary
            %obj.plotLassoCVRestitution_real();
            %obj.plotLassoAMPRestitution();
            
            
        end
        
        
        function getbipolarsignalsandleads(obj)
            %gets signals, leads end electrodenumber from input
            %if falg 'unipolar' is given -> bipolar is calc from unipolar even if
            %bipolar is given
            %if
            if strcmp(obj.signalflag,'unipolar')
                [obj.signal.bi, obj.leads.bi, obj.elecs.bi] = uni2bip(obj.signal.in_uni,obj.leads.in,obj.elecs.in,0);
            elseif strcmp(obj.signalflag,'bipolar')
                obj.signal.bi = obj.signal.in_bi;        %use singal instead of signal.in because a previous filtering can have occured
                obj.leads.bi = obj.leads.in;
                obj.elecs.bi = obj.elecs.in;
            else
                fprintf('Error. signaltype unknown.\n')
            end
            if size(obj.signal.bi,1)<size(obj.signal.bi,2)                     %rows should be time & columns electrodes
                obj.signal.bi = obj.signal.bi';
            end
            if size(obj.leads.bi,1)<size(obj.leads.bi,2)                 %rows should be channels, columns nothing
                obj.leads.bi = obj.leads.bi';
            end
            if size(obj.elecs.bi,1)<size(obj.elecs.bi,2)                 %rows should be channels, columns nothing
                obj.elecs.bi = obj.elecs.bi';
            end
            
        end
        
        
        function splitcatheters(obj)
            [obj.signal.bi_split,obj.leads.bi_split,obj.elecs.bi_split,~] = splitcatheters(obj.signal.bi,obj.leads.bi,obj.elecs.bi,'navx');
            
            if ~isempty(obj.signal.in_uni)
                [obj.signal.uni_split,obj.leads.uni_split,obj.elecs.uni_split,~] = splitcatheters(obj.signal.in_uni,obj.leads.in_uni,obj.elecs.in_uni,'navx');
            end
        end
        
        function getstimchannel(obj)
            obj.stimchan = [];
            obj.stimchan.exclude_chan = [];
            if ~isempty(obj.name)                                                             %if name given try getting from name & from signal
                stimchan = getstimchanfromsig(ECG_Baseline_Removal(obj.signal.bi,obj.samplerate,0.01,0.5),obj.leads.bi,obj.elecs.bi,obj.signalflag); %remove baseline so maximum can be compared
                stimchan2 = getstimchanfromname(obj.name);
                if ~isequal(stimchan.cath,stimchan2.cath)   %if not the same use name, because signals can be tricky
                    fprintf('Stimchan found from signals does not match name of case, check this! Getting stimchan from name\n')
                    obj.stimchan = stimchan2;
                    obj.stimchan.size_chan = size(obj.leads.bi_split.(obj.stimchan.cath),1);
                    obj.stimchan.idpostsplit = find(contains(obj.elecs.bi_split.(obj.stimchan.cath),[num2str(obj.stimchan.el1) '-' num2str(obj.stimchan.el1+1)]));%find(contains(obj.leads.bi_split.SPI,[num2str(obj.stimchan.el1) num2str(obj.stimchan.el1+1)]));
                    %return
                elseif ~isequal(stimchan.el1,stimchan2.el1) && ~isequal(stimchan.el2,stimchan2.el2)
                    fprintf('Stimelectrodes found from signals do not match the ones in name. Using name as default.\n')
                    obj.stimchan = stimchan2;
                    obj.stimchan.size_chan = size(obj.leads.bi_split.(obj.stimchan.cath),1);
                    obj.stimchan.idpostsplit = find(contains(obj.elecs.bi_split.SPI,[num2str(obj.stimchan.el1) '-' num2str(obj.stimchan.el1+1)]));
                    if isempty(obj.stimchan.idpostsplit)
                        obj.stimchan.idpostsplit = str2double(input('Stimechannel is not part of the signals, choose closest matching ID. (Dont count CS if SPI is stimchan)\n','s'));
                    end
                else
                    obj.stimchan = stimchan;
                end
                clearvars stimchan2
            else                                                                                %if name not given, guess stimchan from signal %TODO not robust, needs more love (find bumps after stim in chan next)
                obj.stimchan = getstimchanfromsig(obj.signal,obj.leads.in);
            end
            
            fprintf('Stimchan is: %s\nStimelectrodes are: %.0f and %.0f .\n',obj.stimchan.cath,obj.stimchan.el1,obj.stimchan.el2)
            
            %create non-stimchan as well
            obj.stimchan.non_stimchan = char(setdiff({'CS' 'SPI'},obj.stimchan.cath));
            
            %create excludechannel. For Cs Stim only stimchan, for spiral
            %also left & right neighbour
            obj.stimchan.exclude_chan = [];
            if strcmpi(obj.stimchan.cath,'cs')
                left_of_stimchan = [];
                right_of_stimchan = [];
            else
                if obj.stimchan.idpostsplit-1==0
                    left_of_stimchan = obj.stimchan.size_chan;
                    right_of_stimchan = obj.stimchan.idpostsplit+1;
                elseif obj.stimchan.idpostsplit+1==obj.stimchan.size_chan+1
                    left_of_stimchan = obj.stimchan.idpostsplit;
                    right_of_stimchan = 1;
                else
                    right_of_stimchan = obj.stimchan.idpostsplit + 1;
                    left_of_stimchan = obj.stimchan.idpostsplit - 1;
                end
            end
            obj.stimchan.exclude_chan = unique([obj.stimchan.exclude_chan;obj.stimchan.idpostsplit;right_of_stimchan;left_of_stimchan]);
            
            %find most opposite channel
            possible_chan_ids = [1:obj.stimchan.size_chan];
            non_excluded_ids = possible_chan_ids(~ismember(possible_chan_ids,obj.stimchan.exclude_chan));
            furthest_id = min(diff(sort(non_excluded_ids)));
            
            obj.stimchan.distant_el_id =  furthest_id;
            
        end
        
        function detectEKGpeaks(obj)
            % 1. Pon: Beginning of the P wave
            % 2. Ppeak: Peak of the P wave
            % 3. Poff: End of the P wave
            % 4. QRSon: Beginning of the QRS complex
            % 5. Qpeak: Q peak of the QRS complex
            % 6. Rpeak: R peak of the QRS complex
            % 7. Speak: S peak of the QRS complex
            % 8. QRSoff: End of the QRS complex
            % 9. Lpoint: The middle of the ST segment
            % 10. Pon: Beginning of the T wave
            % 11. Ppeak: Peak of the T wave
            % 12. Poff: End of the T wave
            % 13. Classification of each beat: 1: Normal beat; 2: Ventricular ectopic beat; 3; Supraventricualr ectopic beat; 20: Unclassified
            
            obj.ecg.detection = Annotate_ECG_Multi(obj.ecg.in.ydata(:,1),obj.ecg.in.samplerate);
            
        end
        
        function presegmentsignal(obj)
            f1 = figure;
            figure(f1);
            ax1 = axes;
            plot(ax1,obj.signal.bi_split.(obj.stimchan.cath)(:,obj.stimchan.idpostsplit))
            [x,~] = ginput(2);
            if x(1)<x(2)
                sig_trim_start = floor(x(1));
                sig_trim_end = floor(x(2));
            else
                sig_trim_start = floor(x(2));
                sig_trim_end = floor(x(1));
            end
            
            fieldnam = fieldnames(obj.signal.bi_split);
            for i=1:numel(fieldnam)
                if ~isempty(obj.signal.bi_split.(fieldnam{i}))
                    if sig_trim_end > size(obj.signal.bi_split.(fieldnam{i}),1)
                        sig_trim_end = size(obj.signal.bi_split.(fieldnam{i}),1);
                    end
                    if sig_trim_start < 1
                        sig_trim_start = 1;
                    end
                    obj.signal.bi_split.(fieldnam{i}) = obj.signal.bi_split.(fieldnam{i})(sig_trim_start:sig_trim_end,:);
                end
            end
            
            %unipolar
            if ~isempty(obj.signal.in_uni) %uni_split
                fieldnam = fieldnames(obj.signal.uni_split);
                for i=1:numel(fieldnam)
                    if ~isempty(obj.signal.uni_split.(fieldnam{i}))
                        if sig_trim_end > size(obj.signal.uni_split.(fieldnam{i}),1)
                            sig_trim_end = size(obj.signal.uni_split.(fieldnam{i}),1);
                        end
                        if sig_trim_start < 1
                            sig_trim_start = 1;
                        end
                        obj.signal.uni_split.(fieldnam{i}) = obj.signal.uni_split.(fieldnam{i})(sig_trim_start:sig_trim_end,:);
                    end
                end
            end
            
            %ecg (might have differentsamplerate)
            if ~isequal(size(obj.ecg.in.ydata,1),size(obj.signal.in_bi,1))
                fprintf('Unequal ECG and signal samples\n')
            end
            upsampfactor = obj.ecg.in.samplerate/obj.samplerate;
            if upsampfactor~=1
                for i=1:size(obj.ecg.in.ydata,2)
                    new_vec_dat(:,i) = interp(obj.ecg.in.ydata(:,i),upsampfactor);
                end
                obj.ecg.in.ydata = new_vec_dat;
            end
            obj.ecg.filt.ydata = obj.ecg.in.ydata(sig_trim_start:sig_trim_end,:);
            
            obj.signal.trimid = [sig_trim_start-1 sig_trim_end];
        end
        
        function prefilter_signals(obj)
            %prefilter signals according to what you need here
            %all signals are filtered, so you campare same things later
            
            %======================
            %       Bipolar
            %======================
            obj.signal.bi_split_prefilt = obj.signal.bi_split; %in case no filtering is used take unfilt vals
            
            allcath = fieldnames(obj.signal.bi_split);
            for i=1:numel(allcath)
                sig = obj.signal.bi_split.(allcath{i});
                
                %add your filters here: in form: sig = filter(sig,samplerate,bla)
                sig = ECG_High_Filter(sig,obj.samplerate,obj.init.hp_bi); %obj.init.hp_bi statt 1Hz
                %sig = add_further_filters_her;
                obj.signal.bi_split_prefilt.(allcath{i}) = sig;
            end
            
            %======================
            %       Unipolar
            %======================
            if ~isempty(obj.signal.in_uni)
                obj.signal.uni_split_prefilt = obj.signal.uni_split;
                
                allcath = fieldnames(obj.signal.uni_split);
                for i=1:numel(allcath)
                    sig = obj.signal.uni_split.(allcath{i});
                    
                    %add your filters here: in form: sig = filter(sig,samplerate,bla)
                    sig = ECG_High_Filter(sig,obj.samplerate,1);
                    sig = ECG_Baseline_Removal(sig,obj.samplerate,0.02,0.5); %windowlength without overlap should not be shorter than atrial signal!
                    
                    obj.signal.uni_split_prefilt.(allcath{i}) = sig;
                end
                
                
                %=======================
                % Delete channels at the  shaft of the catheter of 12 and
                % 22 electrode catheters
                %=======================
                %unipolar
                if size(obj.signal.uni_split_prefilt.SPI,2) == 12 || size(obj.signal.uni_split_prefilt,2) == 22
                    obj.signal.uni_split.SPI = obj.signal.uni_split.SPI(:,1:end-2);
                end
            end
            %bipolar
            if size(obj.signal.bi_split_prefilt.SPI,2) == 11 || size(obj.signal.bi_split_prefilt,2) == 21
                obj.signal.bi_split.SPI = obj.signal.bi_split.SPI(:,1:end-2);
            end
        end
        
        
        function getwaveletinfo(obj,k)
            allcath = fieldnames(obj.signal.bi_split_prefilt);
            fieldidofstim = strcmp(obj.stimchan.cath,allcath);
            stimcathname = (allcath(fieldidofstim));
            obj.signal.bi_stimchan = obj.signal.bi_split_prefilt.(stimcathname{1});
            
            stimchan_ids = 1:obj.stimchan.size_chan;
            
            [coeff,score,latent,tsquared,explained,mu] = pca(obj.signal.bi_stimchan(:,stimchan_ids));
            sig = score(:,1);
            
            sig_bi_wavelet = method_wavelet(sig','WaveletShape','bior1.5','samplerate',obj.samplerate,'ThresholdFactor',k,'MinimumInaktivityLength',10,'MinimumAktivityLength',15,'Postprocessing','On'); %smaller k -> more peaks, thresh=k*std 10 or 12
            [~,activseg_wavelet] = getActiveSegmentsFromNLEOInShortWindow(sig_bi_wavelet',obj.samplerate,0.1);   %chose k higher for less (or thinner main) peaks here 
            activseg = repmat(activseg_wavelet,[obj.stimchan.size_chan 1]);
            
            for i=1:obj.stimchan.size_chan
                [channel_peaks_fs{i},~] = getpeaksfromActiveSegmentsNLEO(nleo(obj.signal.bi_stimchan(:,i),obj.samplerate,1,k),activseg{i},obj.samplerate); %active segment peaks to calc distances later
                
                %get most common signallength
                tmp = diff(activseg{i},1,2);         %difference between start and endpoint of active segemnts gives length
                [N,edges] = histcounts(tmp(tmp~=0));
                [~,idx_max_signallength] = max(N);   %max of histogramm
                [~,idx_min_signallength] = min(N);
                mean_signallength_val(i) = mean([edges(idx_max_signallength) edges(idx_max_signallength+1)]);
                lowest_signallength_val(i) = edges(idx_min_signallength+1);
            end
            obj.signal_properties.mean_signallength_fs = floor(mean(mean_signallength_val));              %we will mostly use mean, because assumtion is, that stimuli are large part of signal and contribute mostly to mean
            obj.signal_properties.lowest_signallength_fs = floor(mean(lowest_signallength_val));
            obj.guess_siglength_stim = round(mean([obj.signal_properties.mean_signallength_fs obj.signal_properties.lowest_signallength_fs]));                            %mean siglength is mix of stim & atria & noise, min is probably minimal stim, take mean of the two as guess
            obj.threshold.max_deviation_S1_block_guess = 2*obj.signal_properties.mean_signallength_fs;
            
            obj.signal.activseg = [obj.signal.activseg activseg'];
            obj.signal.channel_peaks_fs = [obj.signal.channel_peaks_fs channel_peaks_fs];
        end
        
        function getBezierinfo(obj,k)
            allcath = fieldnames(obj.signal.bi_split_prefilt);
            fieldidofstim = strcmp(obj.stimchan.cath,allcath);
            stimcathname = (allcath(fieldidofstim));
            obj.signal.bi_stimchan = obj.signal.bi_split_prefilt.(stimcathname{1});
            
            obj.signal.bi_stimchan = obj.signal.bi_split_prefilt.(stimcathname{1});
            
            %Bezier changes method
            stepfcn = zeros(size(obj.signal.bi_stimchan,1),obj.stimchan.size_chan);
            for i=1:obj.stimchan.size_chan
                changes = findchangepts(obj.signal.bi_stimchan(:,i),'Statistic','linear','MinThreshold',0.6);
                diff_changes = diff(changes);
                threshold_ms = 10;
                short_ids = find(diff_changes<threshold_ms*1000/obj.samplerate);
                for ii=1:numel(diff_changes)
                    if diff_changes(ii)<threshold_ms*1000/obj.samplerate
                        stepfcn(changes(ii):changes(ii+1),i) = 1;
                    else
                        
                    end
                end
            end
            [~,activseg] = getActiveSegmentsFromNLEOInShortWindow(stepfcn,obj.samplerate,0.1);   %chose k higher for less peaks here
              
            for i=1:obj.stimchan.size_chan
                [channel_peaks_fs{i},~] = getpeaksfromActiveSegmentsNLEO(nleo(obj.signal.bi_stimchan(:,1),obj.samplerate,1,k),activseg{i},obj.samplerate); %active segment peaks to calc distances later
                
                %get most common signallength
                tmp = diff(activseg{i},1,2);                                      %difference between start and endpoint of active segemnts gives length
                [N,edges] = histcounts(tmp(tmp~=0));
                [~,idx_max_signallength] = max(N);              %max of histogramm
                [~,idx_min_signallength] = min(N);
                mean_signallength_val(i) = mean([edges(idx_max_signallength) edges(idx_max_signallength+1)]);
                lowest_signallength_val(i) = edges(idx_min_signallength+1);
            end
            obj.signal_properties.mean_signallength_fs = floor(mean(mean_signallength_val));              %we will mostly use mean, because assumtion is, that stimuli are large part of signal and contribute mostly to mean
            obj.signal_properties.lowest_signallength_fs = floor(mean(lowest_signallength_val));
            obj.guess_siglength_stim = round(mean([obj.signal_properties.mean_signallength_fs obj.signal_properties.lowest_signallength_fs]));                            %mean siglength is mix of stim & atria & noise, min is probably minimal stim, take mean of the two as guess
            obj.threshold.max_deviation_S1_block_guess = 2*obj.signal_properties.mean_signallength_fs;
            
            obj.signal.activseg = [obj.signal.activseg activseg'];
            obj.signal.channel_peaks_fs = [obj.signal.channel_peaks_fs channel_peaks_fs];
            
        end
        
        function getnleoinfo(obj,k)
            %Uses nleo on all signals, thresholding that dependant on the standard deviation given by k
            %gets a step function and that gives active segments dependant on
            %From there a guess of stimlength in samples is produced and a threshold by how many
            %samples the S1 guess is allowed to differ from the inputtet
            %value.
            
            allcath = fieldnames(obj.signal.bi_split_prefilt);
            fieldidofstim = strcmp(obj.stimchan.cath,allcath);
            stimcathname = (allcath(fieldidofstim));
            obj.signal.bi_stimchan = obj.signal.bi_split_prefilt.(stimcathname{1});
            clearvars stimcathname fieldidofstim allcath
            
            sig = obj.signal.bi_stimchan;
            
            obj.signal.bi_nleo = nleo(sig,obj.samplerate,1,k);                                                  %k is completly useless here
            [~,activseg] = getActiveSegmentsFromNLEOInShortWindow(obj.signal.bi_nleo,obj.samplerate,k);   %chose k higher for less peaks here
           
            chan_id = find(size(activseg)==obj.stimchan.size_chan);
            if chan_id == 1                                                 %columns should always be channels
                activseg = activseg';
            end
            
            for i=1:obj.stimchan.size_chan
                [channel_peaks_fs{i},~] = getpeaksfromActiveSegmentsNLEO(obj.signal.bi_nleo(:,i),activseg{i},obj.samplerate); %active segment peaks to calc distances later
                
                %get most common signallength
                tmp = diff(activseg{i},1,2);                                      %difference between start and endpoint of active segemnts gives length
                [N,edges] = histcounts(tmp(tmp~=0));
                [~,idx_max_signallength] = max(N);              %max of histogramm
                [~,idx_min_signallength] = min(N);
                mean_signallength_val(i) = mean([edges(idx_max_signallength) edges(idx_max_signallength+1)]);
                %highest_signallength_val(i) = edges(idx_max_signallength+1);        %+1 because we want lower edge of next bin to be sure to get the signals in those max bins as well
                lowest_signallength_val(i) = edges(idx_min_signallength+1);
            end
            obj.signal_properties.mean_signallength_fs = floor(mean(mean_signallength_val));              %we will mostly use mean, because assumtion is, that stimuli are large part of signal and contribute mostly to mean
            obj.signal_properties.lowest_signallength_fs = floor(mean(lowest_signallength_val));
            obj.guess_siglength_stim = round(mean([obj.signal_properties.mean_signallength_fs obj.signal_properties.lowest_signallength_fs]));                            %mean siglength is mix of stim & atria & noise, min is probably minimal stim, take mean of the two as guess
            obj.threshold.max_deviation_S1_block_guess = 2*obj.signal_properties.mean_signallength_fs;
            
            obj.signal.activseg = [obj.signal.activseg activseg];
            obj.signal.channel_peaks_fs = [obj.signal.channel_peaks_fs channel_peaks_fs];
        end
        
        function createthreshold_pre_dist(obj)
            %change weighting scheme according to number of channels
            num_of_chan = size(obj.signal.bi_stimchan,2);
            %at least 50% of all nonstim chan haf to be hit and stimchan
            %must have peak!
            obj.threshold.pre_dist = (floor(0.5*(num_of_chan-1))*obj.weight.detected_in_chan+obj.weight.detected_in_stim_chan)+(floor(0.5*((num_of_chan-1)))*obj.weight.hp_detected_in_chan+obj.weight.hp_detected_in_stim_chan);           %2 because we use NLEO & LP Filter, +weight.dist/2 because we want to be above that value but not include the distance weighted probably S1 peaks
            
        end
        
        function createthreshold_dist(obj)
            %change weighting scheme according to number of channels
            num_of_chan = size(obj.signal.bi_stimchan,2);
            %at least 50% of all nonstim chan haf to be hit and stimchan
            %must have peak!
            obj.threshold.dist = (floor(0.5*(num_of_chan-1))*obj.weight.detected_in_chan...
                +obj.weight.detected_in_stim_chan)+(floor(0.5*((num_of_chan-1)))*obj.weight.hp_detected_in_chan+obj.weight.hp_detected_in_stim_chan)...
                +obj.weight.diff...
                +obj.weight.amp-1;           %2 because we use NLEO & LP Filter, +weight.dist/2 because we want to be above that value but not include the distance weighted probably S1 peaks
           % +obj.weight.dist...
                
        end
        
        function clearunneededvars(obj)
            obj.signal.work=[];
            obj.signal.bi=[];
            %obj.signal.bi_split=[];
            obj.leads.bi=[];
        end
        
        function createtimeline(obj)
            %Creates a timeline class wich is basically a cell array
            %it is initiallised with the segments found in the sttimulation
            %channel
            timeline_all = timeline();                          %initialises peaks fromm stimchan
            weightvector = [obj.weight.detected_in_stim_chan,obj.weight.detected_in_chan,obj.weight.hp_detected_in_stim_chan,obj.weight.hp_detected_in_chan,obj.weight.dist];
            thresholdvector = [obj.threshold.deviation_S1_guess];
            timeline_all.initialize_weights(weightvector);      %as of now weights are given as ids referring to this array
            timeline_all.initialize_thresholds(thresholdvector);
            timeline_all.initialize_stimchanid(obj.stimchan.idpostsplit);
            
            for i=1:obj.stimchan.size_chan
                if i == obj.stimchan.idpostsplit
                    timeline_all.addsegments(obj.signal.activseg{i}(:,1),obj.signal.activseg{i}(:,2),obj.signal.channel_peaks_fs{i},1,i,'NLEO');
                else
                    timeline_all.addsegments(obj.signal.activseg{i}(:,1),obj.signal.activseg{i}(:,2),obj.signal.channel_peaks_fs{i},2,i,'NLEO')
                end
            end
            timeline_all.createsegmentids();
            timeline_all.createsupersegments(); % 2 is weightid for detection in another channel, second is threshold id
            timeline_all.createsupersegmentids();
            timeline_all.updateweight();
            obj.timeline = timeline_all();
            
        end
        
        function crosscorrall(obj)
            %crosscorrelating all supersegments with each other in
            %within stimulation channel
            obj.timeline.numofsupersegments;
            cnt_g = 0;
            
            similarlist = cell(0);
            fi_cel=[];
            fj_cel=[];
            
            for i=1:obj.timeline.numofsupersegments
                sig_o{i} = obj.signal.bi_stimchan((obj.timeline.supersegments{i,1}.min_fs:obj.timeline.supersegments{i,1}.max_fs));
                for j=i+1:obj.timeline.numofsupersegments
                    l = obj.timeline.supersegments{i,1}.max_fs-obj.timeline.supersegments{i,1}.min_fs;
                    new_id_val = (obj.timeline.supersegments{j,1}.min_fs:obj.timeline.supersegments{j,1}.min_fs+l);
                    corc = corrcoef(obj.signal.bi_stimchan((obj.timeline.supersegments{i,1}.min_fs:obj.timeline.supersegments{i,1}.max_fs),obj.stimchan.idpostsplit),obj.signal.bi_stimchan(new_id_val,obj.stimchan.idpostsplit));
                    xc_d(i,j) = corc(1,2);
                    
                end
            end
            
            [t_r,t_c] = ind2sub(size(xc_d),find(xc_d>0.9));
            t = [t_r t_c];
            
            for ic = 1:obj.timeline.numofsupersegments
                [r,c] = ind2sub(size(t),find(t==ic));
                g{ic} = unique(t(r,:));
            end
            
            cnt_g = 1;
            while 1
                for ic = 1:size(g,2)
                    logicid(ic) = sum(ismember(g{cnt_g},g{ic}))>0;
                end
                cnt_g = cnt_g+1;
                
            end
            
            
            
            
            for ic = 1:obj.timeline.numofsupersegments
                [r,c] = ind2sub(size(t),find(t==ic));
                g = unique(t(r,:));
                
                if isempty(g)
                    %nogroup{} = g;
                else
                    if isempty(similarlist)
                        similarlist{1} = g;
                    else
                        for icel = 1:size(similarlist,1)
                            fi_cel(icel) = sum(ismember(g,similarlist{icel}))>0;
                        end
                        
                        if sum(fi_cel)>1
                            ids = find(fi_cel==1);
                            comb_g = [];
                            for kk = 1:numel(ids)
                                comb_g = [comb_g; similarlist{ids(kk)}];
                            end
                            comb_g = unique(comb_g);
                            for kk = 1:numel(ids)
                                similarlist{ids(kk)} = [];
                            end
                            
                            similarlist{end+1} = unique([comb_g;g]);
                        elseif sum(fi_cel)==1 %hit detected, add to group
                            similarlist{find(fi_cel==1)} = unique([similarlist{find(fi_cel==1)}; g]);
                            cnt_g =cnt_g+1;
                            fprintf('%d\n',cnt_g)
                        else %create new group
                            similarlist{end+1,1} = g;
                        end
                        
                    end
                end
            end
                        
        end
        
        function high_lp_filt_weighting(obj,k)
            a = ECG_High_Filter(obj.signal.bi_stimchan,obj.samplerate,obj.init.highhp);
            a2 = nleo(a,obj.samplerate,1,k);
            [a2obj.signal,a2obj.activseg] = getActiveSegmentsFromNLEOInShortWindow(a2,obj.samplerate,k);
            for i=1:size(a2obj.signal,2)
                [highfrequ_peaks_fs{i},~] = getpeaksfromActiveSegmentsNLEO(a2,a2obj.activseg{i},obj.samplerate); %active segment peaks to calc distances later
            end
            
            
            for i=1:max(size(highfrequ_peaks_fs))   %channels
                %cnt_j = 0;
                if i==obj.stimchan.idpostsplit      %stimchan
                    cnt_j=0;
                    for j=1:max(size(highfrequ_peaks_fs{i}))
                        %fprintf('i=%d,j=%d\n',i,j)
                        tmp = obj.timeline.checkifpeakfitsexistingsupersegment(highfrequ_peaks_fs{i}(j));
                        if isempty(tmp)
                            %superseg_ids{i}(cnt_j) = 0;
                        elseif sum(size(tmp))<2
                            fprintf('more than one segment hit\n');
                        else
                            cnt_j = cnt_j+1;
                            superseg_ids{i}(cnt_j) = tmp;
                            obj.timeline.addweight2supersegmentbyid(superseg_ids{i}(cnt_j),obj.weight.hp_detected_in_stim_chan)
                        end
                        
                        
                    end
                else                                %otherchan
                    for j=1:max(size(highfrequ_peaks_fs{i}))
                        if isempty(highfrequ_peaks_fs{i})
                            tmp = [];
                        else
                            tmp = obj.timeline.checkifpeakfitsexistingsupersegment(highfrequ_peaks_fs{i}(j));
                        end
                        if isempty(tmp)
                            superseg_ids{i}(j)=NaN;
                        elseif sum(size(tmp))<2
                            fprintf('more than one segment hit\n');
                        else
                            %cnt_j = cnt_j+1;
                            superseg_ids{i}(j)=tmp;
                        end
                        %j=%d\n',i,j)
                        
                        if ~isnan(superseg_ids{i}(j))   %only add weight for existing peaks
                            obj.timeline.addweight2supersegmentbyid(superseg_ids{i}(j),obj.weight.hp_detected_in_chan)
                        end
                    end
                    
                end
            end
        end
        
        function getpeakdistancematrix(obj)
            obj.timeline.update_numofsupersegments;
            if ~isempty(obj.timeline.distance)
                obj.timeline.distance = [];
            end
            
            cnt=0;
            for i = 1:obj.timeline.numofsupersegments
                for j = i+1:size(obj.timeline.supersegments,1)
                    cnt=cnt+1;
                    obj.timeline.distance(cnt,1) = i;
                    obj.timeline.distance(cnt,2) = j;
                    obj.timeline.distance(cnt,3) = abs(obj.timeline.supersegments{i}.peaks(1)-obj.timeline.supersegments{j}.peaks(1));
                    %using peaks, because for always same stimsignal NLEO
                    %peaks shouldnt differ to much -> better than using
                    %mean
                end
            end
        end
        
        function guessS1time(obj)
            
            %use FFT to check if frequency < 300 has been filtered TODO
            %FFT(obj.signal.bi_split.(obj.stimchan.cath),obj.samplerate)
            
            obj.variables.k = obj.init.k - obj.init.step_k; %will be added in while loop again
            itterations_S1_guess = 0;
            obj.signal.channel_peaks_fs = [];
            obj.signal.activseg = [];
            while obj.init.S1_guess_good ~= true
                
                obj.variables.k = obj.variables.k+obj.init.step_k;
                obj.signal.channel_peaks_fs = [];
                obj.signal.activseg = [];
                %TODO allow for plot & manual change & check before running
                %rest of script, since here width of segments is
                %determined
                if strcmpi(obj.variables.method_s1s2detect,'wavelet')
                    obj.getwaveletinfo(0.9);                    %use fixed value for now
                elseif strcmpi(obj.variables.method_s1s2detect,'beziere')
                    obj.getBezierinfo(obj.variables.k);
                else %standard is nleo
                    obj.getnleoinfo(obj.variables.k);             %gets nleo of stimcatheter; higher k -> less peaks, because thresh = k*std(nleosig);
                end
                
                obj.clearunneededvars();
                %try
                    obj.createtimeline();       %as of here supersegments get weighting
                    %obj.timeline.plotsupersegments(obj.signal.bi_stimchan(:,1),0)
                    obj.getpeakdistancematrix
                    
                    for l = 1:obj.threshold.max_neighbour
                        logic_row_id = diff([obj.timeline.distance(:,1),obj.timeline.distance(:,2)],1,2)==l;
                        neighbour_dist_fs = obj.timeline.distance(logic_row_id,3);
                        
                        %get maximum of histogramm to guess S1 distance
                        samplerate_round = round(obj.samplerate);
                        [N,edges] = histcounts(neighbour_dist_fs,samplerate_round);
                        %[N,edges] = histcounts(peaks_dist(peaks_dist~=0),2*obj.samplerate);
                        [val(l),idx] = max(N);
                        multiple_max_val = find(N==val(l));
                        if size(multiple_max_val,2) > 1
                            idx = multiple_max_val(end);                %take last one as best guess
                        end
                        bin_cent(l) = mean([edges(idx) edges(idx+1)]);
                        %old code
                        %S1_guess(l,1) = bin_cent(l)/l;
                    end
                    %%
                    cnts1=0;
                    for ir=1:size(bin_cent,2)
                        for ic=ir+1:size(bin_cent,2)
                            [ismult(ir,ic),multfac(ir,ic)] = ismultiple(bin_cent(ir),bin_cent(ic),obj.threshold.deviation_S1_guess);
                            if ismult(ir,ic)==1
                                cnts1=cnts1+1;
                                S1_guess(cnts1) = max([bin_cent(ir),bin_cent(ic)]/multfac(ir,ic));
                            end
                        end
                    end
                    S1_check_sum = sum(ismult(:)); %subtr diagonal; always 1)
%                     for l = 1:length(bin_cent)
%                         remain(l,:) = mod(bin_cent,bin_cent(l));
%                     end
%                     for ir=1:size(remain,1)
%                         for ic=1:size(remain,2)
%                             if ic==ir
%                                 remain(ir,ic) = nan;
%                             end
%                         end
%                     end
%                     figure
%                     histogram(remain,samplerate_round)

%                     %check 1: check neighbouring harmonics 
%                     for i=1:obj.threshold.max_neighbour
%                         for j=1:obj.threshold.max_neighbour
%                             S1_check(i,j) = abs(S1_guess(i)-S1_guess(j));
%                             if i>j
%                                 S1_check(i,j)=0;
%                             end
%                         end
%                     end
%                     S1_check2 = S1_check(S1_check~=0);
%                     S1_check_sum = S1_check2<=obj.threshold.deviation_S1_guess;%/100*mean(S1_guess);
                    
                    %check 2: check mean stimlength:
                    if obj.timeline.getmediansiglength - obj.threshold.stimuluslength_fs <= obj.threshold.stimuluslength_fs_deviation
                        fprintf('mean stimulussignallength: %f ms',obj.timeline.getmeansiglength/obj.samplerate*1000)
                        obj.init.S1_guess_good=true;
                        S1_guess_final_exact_fs = obj.timeline.getmediansiglength;
                    end
                    
                    if sum(S1_check_sum)>floor(obj.threshold.max_neighbour/2)   %2 because we use 50% as right detection
                        %old code
                        %[id1 id2] = ind2sub([obj.threshold.max_neighbour obj.threshold.max_neighbour],find(min(S1_check2)==S1_check)); %get the two ids closest together
                        %S1_guess_final_exact_fs = S1_guess(id1(1));                                                 %for now just take any found value
                        [v,c]=kmeans([1:length(S1_guess); S1_guess]',2);
                        most_common = mode(v);
                        %S1_guess_final_exact_fs = c(most_common,2); 
                        S1_guess_final_exact_fs = S1_guess(find(v==most_common,1,'first')); 
                        
                        obj.init.S1_guess_good=true;
                        fprintf('S1 guess found: %f ms.\n',S1_guess_final_exact_fs/obj.samplerate*1000)%S1_guess_final_exact_fs/obj.samplerate
                        obj.S1_time = S1_guess_final_exact_fs;
                    else
                        obj.init.S1_guess_good=false;
                        itterations_S1_guess = itterations_S1_guess + 1;
                        fprintf('Trying new itteration of weighting.\n')
                        %will jump to beginning of while loop and increase k, allowing less
                        %peaks and hopefully reducing the error
                    end
                    
                    if itterations_S1_guess == obj.threshold.max_itterations_S1_guess
                        fprintf('No good S1 guess could be obtained.\n')
                        return
                    end
                    
                %catch
                %end
            end %Trying to guess S1 while loop
            
            %refining S1_variance threshold if new estimate is smaller
            if 0.05*obj.S1_time < obj.threshold.deviation_S1_guess %TODO change this later by itterating through finer binning
                obj.threshold.deviation_S1_guess = 0.05*obj.S1_time;
            end
            %obj.timeline.plotsupersegmentweights(obj.signal.bi_stimchan)
            
            
            %compare against notated values
            if ~isempty(obj.classcheck)
                fprintf('Expected S1 is %d ms. \n',obj.classcheck.s1time)
                variance_S1_time =round(abs(S1_guess_final_exact_fs/obj.samplerate*1000-obj.classcheck.s1time));
                if variance_S1_time > 25 %25 is just a guess, set threshold however you deem fit
                    fprintf('Difference in S1 is %d ms. This value will be set manually to %d .\n',variance_S1_time,obj.classcheck.s1time)
                    fprintf('This could indicate bad signal quality. \n')
                    obj.S1_time = obj.classcheck.s1time/1000*obj.samplerate;
                end
            end
            
%             %===PLOT for paper===
%              f_temp = figure;
%              histogram(obj.timeline.distance(:,3),samplerate_round);
%             xlim([0 5000]);
%             ylim([0 100])
%             save_name = 'f_hist_noise';
%             savefig(f_temp,['Pics/' save_name '.fig'])
%             print(f_temp,['Pics/' save_name],'-dpng','-r400')
            %
        end
        
        function getpointsinbetweenblocks(obj)
              manual_seg = segmentegmtool(sig_with_buffer,obj.samplerate);
              segments_start = floor(manual_seg(1:end-1,1)) - numel(buffer_left) +1;
              segments_end = floor(manual_seg(2:end,1)) - numel(buffer_left) -1;
              
              %take one block
              for i=1:size(segments_start)
                  cutblock{i} = sig_with_buffer(segments_start(i):segments_end(i));
                  [sigmoid{i},~]=getActiveSegmentsFromNLEOInShortWindow(nleo(cutblock{i},obj.samplerate,1),obj.samplerate);
                  [S1_block_guess_check(i),~,~] = checkguessS1blockcorr(obj.S1_time,6,obj.threshold.s1_block_guess_variance,obj.guess_siglength_stim,sigmoid{i});
              end
              guess_s1_in_block = mode(S1_block_guess_check);
              fprintf('Initial guess of blocknum: %i \n',guess_s1_in_block)
              %findpeaks with diff + saks
              
              
              %sigmoid_active = zeros(numel(sig_with_buffer),1);
              %for i=1:numel(segments_start)
              %    sigmoid_active(segments_start(i):segments_end(i))=1;
              %end
              
              
              [S1_block_guess_check,S1_block_guess_num,loc_of_peaks_fs] = checkguessS1blockcorr(obj.S1_time,6,obj.threshold.s1_block_guess_variance,obj.guess_siglength_stim,sigmoid{1});
         end
        
        function guessblocklocations(obj)
            
            [~,score,~,~,~,~] = pca(obj.signal.bi_stimchan);
            %buffer_left = zeros(5000,1);
            %buffer_right = zeros(5000,1);
            sig = score(:,1);
            %sig_with_buffer = [buffer_left;score(:,1);buffer_right]; %add zeros to avoid incomplete peaks
            
            sig_squared = abs(sig).^2;
            diff_sig = diff(sig_squared); %Squared because then more difference between baseline & peaks
            
            k = 5;
            peak_height_thresh = median(diff_sig) + k * std(diff_sig);
            min_distance_fs = 0.1*obj.S1_time;% min(obj.classcheck.s2times); %90% of s1time
            locs = [];
            for i=size(diff_sig,2)
                [~,locs{i,1}] = findpeaks(diff_sig,'MinPeakHeight',peak_height_thresh,'MinPeakDistance',min_distance_fs);
            end
            
            sigmoidsig = zeros(numel(diff_sig),1);
            for i=1:size(locs,1)
                [cluster_id,~] = kmeans(diff(locs{i}),2,'Start',[min(diff(locs{i})); max(diff(locs{i}))]);
                small_dist_cluster_id = mode(cluster_id);
                large_dist_cluster_id = mode(cluster_id(cluster_id~=small_dist_cluster_id));
                
                small_dist_in_loc = find(cluster_id==small_dist_cluster_id);
                large_dist_in_loc = find(cluster_id==large_dist_cluster_id);
                
                S1_in_block = mode(diff(large_dist_in_loc)) - 1; %-1 to substract S2 stim
                %Quickfix if fits few segments are missdetected - DOESNT
                %WORK
                while large_dist_in_loc(1)-(S1_in_block+1)>0
                    fprintf('This should never happen\n')
                    large_dist_in_loc = large_dist_in_loc(1:end-1) + 1;
                    large_dist_in_loc = [large_dist_in_loc(1)-(S1_in_block+1); large_dist_in_loc];
                    small_dist_in_loc(find(small_dist_in_loc==large_dist_in_loc(1) ))= [];
                end
                
                num_of_blocks = numel(large_dist_in_loc) + 1;
                for j=1:num_of_blocks %start
                    if j==1 && cluster_id(1)==small_dist_cluster_id
                        segments_start(j,1) = locs{i}(1);
                        segments_end(j,1) = locs{i}(large_dist_in_loc(1));
                        sigmoidsig(segments_start(j,1):segments_end(j,1)) = 1;
                        mean_bloc_width_fs(j,1) = segments_end(j,1) - segments_start(j,1);
                    elseif j==1 && cluster_id(1)==large_dist_cluster_id
                        %find next small cluster
                        first_small = find(small_dist_in_loc>large_dist_in_loc(j),1,'first');
                        segments_start(j,1) = locs{i}(small_dist_in_loc(first_small));
                        segments_end(j,1) = locs{i}(large_dist_in_loc(j+1));
                        sigmoidsig(segments_start:segments_end(j,1)) = 1;
                        mean_bloc_width_fs(j,1) = segments_end(j,1) - segments_start(j,1);
                    elseif j==num_of_blocks %end
                        if locs{i}(large_dist_in_loc(j-1)+1)+ceil(median(mean_bloc_width_fs)) <= numel(diff_sig)
                            segments_start(j,1) = locs{i}(large_dist_in_loc(j-1)+1);
                            segments_end(j,1) = locs{i}(large_dist_in_loc(j-1)+1)+ceil(median(mean_bloc_width_fs));
                            sigmoidsig( segments_start(j,1):segments_end(j,1) ) = 1;
                        else
                            segments_start(j,1) = locs{i}(large_dist_in_loc(j-1)+1);
                            segments_end(j,1) = numel(diff_sig)-1; %-1 so sigmoid will always end on zero
                            sigmoidsig( segments_start(j,1):segments_end(j,1) ) = 1;
                            mean_bloc_width_fs(j,1) = mean(mean_bloc_width_fs);
                        end
                    else
                        segments_start(j,1) = locs{i}(large_dist_in_loc(j-1)+1);
                        segments_end(j,1) = locs{i}(large_dist_in_loc(j));%locs{i}(large_dist_in_loc(j))+1
                        sigmoidsig(segments_start(j,1):segments_end(j,1) ) = 1;
                        mean_bloc_width_fs(j,1) = segments_end(j,1) - segments_start(j,1);
                    end
                    
                end
            end
            
            buff_perc = 0.05; %buffer in percent since sigmoid found by combining peaks is very close to actual stimulussignals
            buffer_fs = floor(buff_perc*mean(mean_bloc_width_fs));
            
            segments_start = segments_start - buffer_fs;
            segments_end = segments_end + buffer_fs;
            segments_start(segments_start<1) = 1;
            segments_end(segments_end>numel(sig)) = numel(sig);
            
            for ii=1:length(segments_start)
                sigmoidsig(segments_start(ii):segments_end(ii)) = 1;
            end
            
            fprintf('%i blocks were detected. \n',num_of_blocks)
             
            %             segstart = find(diff(sigmoid)==1);
            %             segend = find(diff(sigmoid)==-1);
%             if length(segstart)>length(segend)
%                 segend = find(diff([sigmoid; 0])==-1);
%             end
%             if segend(1)<segstart(1)
%                 segstart=[1; segstart];
%             end
            
           
            if obj.settings.manualmode
            if obj.settings.ploton
                figure
                hold on
                plot(sig)
                plot(sigmoidsig)
            end
                answer = input('Do you wish to edit this manually? [y,n]: ','s');
            else
                answer = 'n';
            end
            if isempty(answer) answer='n'; end
            if strcmp(answer,'n')
                obj.S1_blocknum = num_of_blocks;
                obj.S1_blocknum_locs = [segments_start,segments_end,mean([segments_start,segments_end],2)];
                obj.numof_S1_in_one_block = S1_in_block;
            else
                manual_seg = segmentegmtool(sig,obj.samplerate);
                segments_start = floor(manual_seg(1:2:end,1));
                segments_end = floor(manual_seg(2:2:end,1));
                obj.numof_S1_in_one_block = input('How many S1 are there in one block? [integer num]: ');
                obj.S1_blocknum = numel(segments_start);
                obj.S1_blocknum_locs = [segments_start,segments_end,mean([segments_start,segments_end],2)];
            end
        end
        
        function guessnumofstimuliinblock(obj)
            blocklength_fs = obj.S1_blocknum_locs(:,2)-obj.S1_blocknum_locs(:,1);
            stimnum_guess = floor(blocklength_fs / obj.S1_time);
            
            obj.numof_stimuli_in_one_block = mode(stimnum_guess);
            obj.numof_S1_in_one_block = mode(stimnum_guess) - 1; %since we know we only get 1 S2 beat
            fprintf('Found amount of S1 peaks in blocks: %i .\n',obj.numof_S1_in_one_block)
        end
        
        function deletesegmentsoutsideblocks(obj)
            cnt = 0;
            idlist = [];
            for i=1:obj.S1_blocknum
                for j=1:obj.timeline.numofsupersegments
                    if obj.timeline.supersegments{j}.peaks(1) > obj.S1_blocknum_locs(i,1) && obj.timeline.supersegments{j}.peaks(1) < obj.S1_blocknum_locs(i,2)
                        cnt=cnt+1;
                        idlist(cnt) = obj.timeline.supersegments{j}.segment_id;
                    end
                end
            end

            for i=1:obj.timeline.numofsupersegments
                if ~ismember(obj.timeline.supersegments{i}.segment_id,idlist)
                    obj.timeline.supersegments{i}.weight = 0;
                end
            end
        end
        
        function guessnumofS1inblock(obj)
            
            all_weights = obj.timeline.getuniqueweights;
            obj.variables.deviation_S1_block_guess = obj.threshold.deviation_S1_block_guess;
            
            if size(all_weights,2)==2 %if only 2 weights, all possible peaks were found
                all_weights = max(all_weights);
                good_supersegids = obj.timeline.getsupersegids_larger_than_weight(all_weights-1)';
                reduced_distance = obj.timeline.createreduceddistmat(good_supersegids);
                rows = find(abs(diff(reduced_distance(:,1:2),1,2)==1));
                reduced_distance2 = reduced_distance(rows,:);
                [kmeansres_idx kmeansres_center] = kmeans(reduced_distance2(:,3),3,'Start',[min(reduced_distance2(:,3));1/3*(max(reduced_distance2(:,3))-min(reduced_distance2(:,3)))+min(reduced_distance2(:,3)) ;max(reduced_distance2(:,3))]);
                %find out wich one of the clusters is s1 (maximum)
                counts(1)=sum(kmeansres_idx==1);
                counts(2)=sum(kmeansres_idx==2);
                counts(3)=sum(kmeansres_idx==3);
                s1_count_id = find(counts==max(counts));
                counts1=0;
                max_cnt_save=0;
                for i=size(kmeansres_idx,1):-1:1 %reverse loop cause lower S2 better distinguishable
                    if kmeansres_idx(i)==s1_count_id
                        counts1 = counts1+1;
                        s1_id_list(counts1) = i;
                    else
                        max_cnt_save(end+1)=i; %These are S2 ids
                    end
                end
                changes2s2 = abs(diff(max_cnt_save(2:end)));
                most_common_num = mode(changes2s2); %this diff counts all spaces between this and next s2 peak so no need to add 1
                guess_numof_S1_in_block = most_common_num;
                if mod(size(good_supersegids,1),guess_numof_S1_in_block+1)~=0
                    fprintf('Shit.\n')
                    S1_block_guess_num = sum(changes2s2==most_common_num)+2;
                else
                    S1_block_guess_num = size(good_supersegids,1)/(guess_numof_S1_in_block+1);
                end
                
                %count peaks to find s2 and locations
                count=0;
                s1_id_list_final=0;
                s2_id_list_final=0;
                for i=1:size(good_supersegids,1)
                    %count = count + 1;
                    if mod(i,guess_numof_S1_in_block+1)==0
                        s2_id_list_final(end+1) = good_supersegids(i); %my head hurts so probably this whole code sucks
                    else
                        s1_id_list_final(end+1) = good_supersegids(i);
                    end
                end
                s1_id_list_final = s1_id_list_final(2:end);
                s2_id_list_final = s2_id_list_final(2:end);
                
                
                obj.init.S1_block_guess_good = true;
                obj.numof_S1_in_one_block = guess_numof_S1_in_block;
                obj.S1_blocknum = S1_block_guess_num;
                %obj.S1_blocknum_locs = loc_of_peak_fs;
                
                all_peak_ids_unique = s1_id_list_final;
                
                for i=1:max(size(s1_id_list_final))
                    obj.timeline.supersegments{s1_id_list_final(i)}.addsegtype('S1');
                end
                for i=1:max(size(s2_id_list_final))
                    obj.timeline.supersegments{s2_id_list_final(i)}.addsegtype('S2');
                end
                
            else
                obj.createthreshold_pre_dist;
                all_weights = all_weights(all_weights>obj.threshold.pre_dist);
                weight_cnt = 0;
                
                if isempty(all_weights)
                    fprintf('No saatisfactory S1 Block found. Exiting...\n');
                    return
                end
                
                while obj.init.S1_block_guess_good~=true
                    %with each loop increase weighting till max_weight, then
                    %increase threshold for S1 signallength deviation
                    
                    weight_cnt = weight_cnt + 1;
                    
                    if weight_cnt>=max(size(all_weights))
                        fprintf('Maximum weight reached. Increasing threshold & retrying.\n')
                        obj.variables.deviation_S1_block_guess = obj.variables.deviation_S1_block_guess + obj.threshold.increase_S1_block_deviation; %increases deviation for stimulussignallength
                        weight_cnt = 1;               %reset weightcounter
                        if obj.variables.deviation_S1_block_guess == obj.threshold.max_deviation_S1_block_guess
                            fprintf('No saatisfactory S1 Block found. Exiting...\n');
                            return
                        end
                    end
                    weight_S1_block_guess = all_weights(weight_cnt);
                    
                    void_supersegids = obj.timeline.getsupersegids_smaller_than_weight(weight_S1_block_guess)';
                    good_supersegids = obj.timeline.getsupersegids_larger_than_weight(weight_S1_block_guess-1)';
                    
                    
                    
                    %tmp = 1:obj.timeline.numofsupersegments;
                    %good_supersegids = tmp(~ismember(tmp,void_supersegids));
                    
                    if isempty(void_supersegids)
                        reduced_distance = obj.timeline.distance();
                    else
                        %reduced_distance = obj.timeline.distance(good_supersegids,:);
                        reduced_distance = obj.timeline.createreduceddistmat(good_supersegids);
                    end
                    
                    %find ids that fit distance
                    [s1_row,~] = find(reduced_distance(:,3)<obj.S1_time+obj.variables.deviation_S1_block_guess & reduced_distance(:,3)>obj.S1_time-obj.variables.deviation_S1_block_guess & reduced_distance(:,3)~=0);
                    
                    %only works if many peaks inbetween
                    if ~isempty(s1_row)
                        all_peak_ids = reduced_distance(s1_row,1:2);
                        [all_peak_ids_unique,~,ic] = unique(all_peak_ids);
                        
                        %if ids are shared
                        N = numel(all_peak_ids_unique);
                        count = zeros(N,1);
                        for k = 1:N
                            count(k) = sum(all_peak_ids(:)==all_peak_ids_unique(k));
                        end
                        %disp([ C(:) count ]);
                        max_count = max(count);
                        
                        %find longest row of max_count
                        cnt = 0;
                        max_row_cnt = 0;
                        for i=1:size(count,1)
                            
                            if count(i)==max_count
                                cnt = cnt + 1;
                            else
                                max_row_cnt(end+1) = cnt;
                                cnt = 0;
                            end
                        end
                        
                        
                        
                        
                        %TODO not a good method
                        if max_row_cnt == 0 %all ids are exactly one count apart
                            %Check how often multiple distances, when count goes
                            %down is bad
                            for i=1:100
                                [s1_row,~] = find(reduced_distance(:,3)<i*(obj.S1_time+obj.variables.deviation_S1_block_guess) & reduced_distance(:,3)>i*(obj.S1_time-obj.variables.deviation_S1_block_guess) & reduced_distance(:,3)~=0);
                                counter(i) = size(s1_row,1);
                                ismin = islocalmin(counter);
                                if sum(ismin)==1
                                    break
                                end
                            end
                            if sum(ismin)==1
                                guess_numof_S1_in_block = find(ismin==1);
                                most_frequent_num = guess_numof_S1_in_block; %to prevent failing in if statement later
                            else
                                guess_numof_S1_in_block = NaN;
                            end
                            %figure
                            %plot(counter)
                        else
                            most_frequent_num = mode(max_row_cnt(max_row_cnt~=0)) + 2; %TODO Think about this long and hard
                            guess_numof_S1_in_block = max(max_row_cnt) + 2; %+2 because diff and last peak
                        end
                        
                        
                    else
                        guess_numof_S1_in_block = NaN; %+2 because diff and last peak
                    end
                    
                    %check S1 block guess with crosscorr
                    if ~isempty(guess_numof_S1_in_block) && guess_numof_S1_in_block > 1 && ~isnan(guess_numof_S1_in_block)
                        %check with crosscorr
                        block_check_signal = obj.signal.step(:,obj.stimchan.idpostsplit);                    %TODO for now only check with stimuluschannel, maybe later with all chans
                        if obj.guess_siglength_stim <= 60 %TODO for now
                            obj.guess_siglength_stim = 100;
                        end
                        %if abs(obj.S1_time-1222) > 400 %TODO for now
                        %    obj.S1_time = 1222;
                        %end
                        [S1_block_guess_check_crosscor,S1_block_guess_num,loc_of_peak_fs] = checkguessS1blockcorr(obj.S1_time,guess_numof_S1_in_block,obj.threshold.s1_block_guess_variance,obj.guess_siglength_stim,block_check_signal);
                        
                        %                         figure
                        %                         hold on
                        %                         plot(obj.signal.bi_stimchan(:,obj.stimchan.idpostsplit))
                        %                         plot(loc_of_peak_fs,repmat(1,numel(loc_of_peak_fs)),'r*','LineWidth',2)
                        
                        if guess_numof_S1_in_block==S1_block_guess_check_crosscor || most_frequent_num==S1_block_guess_check_crosscor %|| most_frequent_num==guess_numof_S1_in_block
                            fprintf('Num of S1 in a Block guess was %.f, check found %.f, similarity: %s\n', guess_numof_S1_in_block ,S1_block_guess_check_crosscor,mat2str(guess_numof_S1_in_block==S1_block_guess_check_crosscor) )
                            %guess_numof_S1_in_block = true;
                            obj.init.S1_block_guess_good = true;
                            obj.numof_S1_in_one_block = guess_numof_S1_in_block;
                            obj.S1_blocknum = S1_block_guess_num;
                            block_radius_fs = floor((obj.numof_S1_in_one_block-1)/2 * obj.S1_time); %block_radius
                            obj.S1_blocknum_locs = [loc_of_peak_fs-block_radius_fs loc_of_peak_fs+block_radius_fs loc_of_peak_fs]; %peaks only in 3rd dimension
                            obj.S1_blocknum_locs(obj.S1_blocknum_locs<0)=0;
                        else
                            fprintf('Num of S1 in a Block guess was %.f, check found %.f, similarity: %s\n', guess_numof_S1_in_block ,S1_block_guess_check_crosscor,mat2str(guess_numof_S1_in_block==S1_block_guess_check_crosscor) );
                            %                             answer = questdlg('Do you want to accept this anyway?','Yes','No');
                            %                             switch answer
                            %                                 case 'Yes'
                            %                                     obj.init.S1_block_guess_good = true;
                            %                                     obj.numof_S1_in_one_block = guess_numof_S1_in_block;
                            %                                     obj.S1_blocknum = S1_block_guess_num;
                            %                                     obj.S1_blocknum_locs = loc_of_peak_fs;
                            %                                 case 'No'
                            %                                     %return
                            %                            end
                        end
                    end
                end %while S1_block_guess_good ~= true
                
            end
        end %function end
        
        function adddistanceweighting(obj)
            %obj.createthreshold_pre_dist;
            [coeff,score,latent,tsquared,explained,mu] = pca(obj.signal.bi_stimchan);
            %buffer_left = zeros(5000,1);
            %buffer_right = zeros(5000,1);
            sig = score(:,1);
            %sig_with_buffer = [buffer_left;score(:,1);buffer_right]; %add zeros to avoid incomplete peaks
            
            sig_squared = abs(sig).^2;
            diff_sig = diff(sig_squared); %Squared because then more difference between baseline & peaks
            diff_diff_sig = diff(diff_sig);
            %dddd_sig = diff(diff_diff_sig);
            
            k = 5;
            peak_height_thresh = median(diff_sig) + k * std(diff_sig);
            min_distance_fs = 0.1*obj.S1_time;% min(obj.classcheck.s2times); %90% of s1time
            for i=size(diff_sig,2)
                [~,locs{i,1}] = findpeaks(diff_sig,'MinPeakHeight',peak_height_thresh,'MinPeakDistance',min_distance_fs);
            end
            
            for i=1:obj.timeline.numofsupersegments
                for j=1:numel(locs{1,1})
                    if obj.timeline.supersegments{i,1}.min_fs<locs{1,1}(j) && obj.timeline.supersegments{i,1}.max_fs>locs{1,1}(j)
                        obj.timeline.supersegments{i,1}.weight = obj.timeline.supersegments{i,1}.weight + obj.weight.dist;
                    end
                end
            end
            obj.timeline.plotsupersegmentweights(obj.signal.bi_stimchan)
            
        end
        
        function addthresholdweighting(obj)
            %add weighting based on voltage threshold
            %obj.createthreshold_pre_dist;
            [coeff,score,latent,tsquared,explained,mu] = pca(obj.signal.bi_stimchan);
            sig = score(:,1);
            
            %manual_seg = segmentegmtool(sig,obj.samplerate);
            %peak_height_thresh = floor(manual_seg(1:1:end,1));
            %histogramm automatic detection using:
            %findpeaks(hist(sig.^2,100)) %TODO
            peak_height_thresh = max(sig)*0.9;%1.99;
            
            min_distance_fs = 0.1*obj.S1_time;% min(obj.classcheck.s2times); %90% of s1time
            [~,locs{1,1}] = findpeaks(sig,'MinPeakHeight',peak_height_thresh,'MinPeakDistance',min_distance_fs);
            
            
            for i=1:numel(locs{1,1})
                for j=1:obj.timeline.numofsupersegments
                    diff_tmp(j,i) = abs(obj.timeline.supersegments{j}.peaks(1)-locs{1,1}(i));
                end
            end
            [~,minidx] = min(diff_tmp,[],1); %saves index of supersegment
            
            for i=1:numel(minidx)
                obj.timeline.supersegments{minidx(i),1}.weight = obj.timeline.supersegments{minidx(i),1}.weight + obj.weight.amp;
            end
%             for i=1:obj.timeline.numofsupersegments
%                 for j=1:numel(locs{1,1})
%                     if obj.timeline.supersegments{i,1}.min_fs<locs{1,1}(j) && obj.timeline.supersegments{i,1}.max_fs>locs{1,1}(j)
%                         obj.timeline.supersegments{i,1}.weight = obj.timeline.supersegments{i,1}.weight + 1000;
%                     end
%                 end
%             end
            if obj.settings.ploton
                obj.timeline.plotsupersegmentweights(obj.signal.bi_stimchan);
            end
            
        end
        
        function adddistanceweighting_old(obj)
            %takes peaks (preselected by weight) and checks theire distance
            %if the distance to the last peak is S1distance and to the next
            %is large then it should be an S2 peak
            %differentiation between large and not large is done via kmeans
            obj.createthreshold_pre_dist;
            all_peak_ids_unique = obj.timeline.getsupersegids_larger_than_weight(obj.threshold.pre_dist);
            
            fprintf('Expecting %i segments. \n',obj.S1_blocknum*(obj.numof_S1_in_one_block+1) ) %+1 S2
            reduced_distance = obj.timeline.createreduceddistmat(all_peak_ids_unique);
            
            mapnew2oldid = [ all_peak_ids_unique;[1:1:max(size(all_peak_ids_unique))] ]'; %col1: supersegid - col2: enum
            change_id = diff(mapnew2oldid,1,2)~=0;
            change_list = mapnew2oldid(change_id,:);
            for i=1:size(change_list,1)
                reduced_distance(reduced_distance(:,1:2)==change_list(i,1)) = change_list(i,2);
            end
            
            %rough reduction
            %find neighbouring ids with dist<2*s1
            neighb_ids = diff([reduced_distance(:,1) reduced_distance(:,2)],1,2)==1;
            neighb_distances = reduced_distance(neighb_ids,:);
            
            %create prelimlist with S1=1 Method1:
            %forward_dist_id = find(neighb_distances(:,3)<=obj.S1_time+obj.threshold.deviation_S1_block_guess & neighb_distances(:,3)>=obj.S1_time-obj.threshold.deviation_S1_block_guess);
            %prelim_list_s1 = unique(neighb_distances(forward_dist_id,1:2));
            
            %create prelimlist with S1=1 Method2:
            numofclust = 1:2; %We want to distinguish the peaks and rest
            [cluster_list, ~] = kmeans(neighb_distances(:,3),numofclust(end),'Start',[min(neighb_distances(:,3));max(neighb_distances(:,3))]);
            min_cluster = mode(cluster_list); %this should be min value
            max_cluster = setdiff(numofclust,min_cluster);
            
            un = unique(neighb_distances);
            ids_blocks = un(cluster_list==min_cluster);
            %ids_block_end = neighb_distances(find(cluster_list==max_cluster),3); works nicely!
            prelim_list_s1 = unique(ids_blocks);
            
            %converting to logical list
            logical_list = zeros(max(size(all_peak_ids_unique)),1);
            logical_list(prelim_list_s1) = 1;
            expected_sequence_one_block = [ones(obj.numof_S1_in_one_block,1)' 0];
            s1_blocks_start(:,2) = strfind(logical_list',expected_sequence_one_block);
            
            %transition_ids = find(diff(forward_dist_id)>1);
            %bad_dist_ids = prelim_list_s1(transition_ids+1);
            
            %obj.timeline.plotsupersegmentweights(obj.signal.bi_split.(obj.stimchan.cath))
            
            if size(s1_blocks_start,2) >= obj.S1_blocknum
                %get centers
                conv2oldid = all_peak_ids_unique(s1_blocks_start);
                for i=1:obj.S1_blocknum
                    %Method 1
                    if mod(obj.numof_S1_in_one_block,2)==0 %if even
                        s1_blocks_starts(i,1) = round(mean(obj.timeline.supersegments{conv2oldid(i)}.peaks));
                        s1_blocks_centers(i,1) = round(mean([obj.timeline.supersegments{conv2oldid(i)+obj.numof_S1_in_one_block/2-1}.peaks,obj.timeline.supersegments{conv2oldid(i)+obj.numof_S1_in_one_block/2}.peaks]));
                        s1_blocks_ends(i,1) = round(mean(obj.timeline.supersegments{conv2oldid(i)+(obj.numof_S1_in_one_block-1)}.peaks));
                    else
                        fprintf('Check me for errors\n')
                        s1_blocks_starts(i,1) = round(mean(obj.timeline.supersegments{conv2oldid(i)}.peaks));
                        s1_blocks_centers(i,1) = round(mean(obj.timeline.supersegments{conv2oldid(i)+(obj.numof_S1_in_one_block-1)/2}.peaks));
                        s1_blocks_ends(i,1) = round(mean(obj.timeline.supersegments{conv2oldid(i)+(obj.numof_S1_in_one_block-1)}.peaks));
                    end
                    %Method 2 (probably not so good)
                    %s1_block_centers(i,1) = round(mean(obj.timeline.supersegments{conv2oldid(i)}.peaks) + (obj.numof_S1_in_one_block-1)/2 * obj.S1_time);
                    %s1_block_ends(i,1) = round(mean(obj.timeline.supersegments{conv2oldid(i)}.peaks) + (obj.numof_S1_in_one_block-1) * obj.S1_time);
                end
                obj.S1_blocknum_locs = [s1_blocks_starts,s1_blocks_ends,s1_blocks_centers];
                
                %dismiss centers outsiede range
            elseif isempty(obj.S1_blocknum_locs) || size(obj.S1_blocknum_locs,1)~=obj.S1_blocknum
                fprintf('strfind went wrong.\n')
                return
            end
            %expected_sequence_one_block = [ones(obj.numof_S1_in_one_block,1)' 2 ];
            %expected_sequence_two_blocks = [zeros(obj.numof_S1_in_one_block,1)' 1 zeros(obj.numof_S1_in_one_block,1)'];
            %expected_sequence_full = [zeros(obj.numof_S1_in_one_block,1)' 1 zeros(obj.numof_S1_in_one_block,1)'];
            
            
            %get S1 blocks
            %
            %
            %             cnts1 = 0;
            %             for i=1:size(prelim_list_s1,1)
            %                 cnts1=cnts1+1;
            %                 if cnts1==obj.numof_S1_in_one_block
            %                     cnts1=0;
            %                     blocks{i,1}=1;
            %                 end
            %             end
            
            %         transition_ids = find(diff(prelim_list_s1)>1);
            
            
            %             if isempty(obj.S1_blocknum_locs) %if for some reason this should be empty
            %
            %             end
            %             not_s1 = neighb_distances(~ismember(all_peak_ids_unique,prelim_list_s1));
            %             prelim_list_s2 = 1;
            %             %split into prelim blocks
            %
            %             %fine reduction
            %             kmeans(reduced_distance(:,3),4);
            %             expected_num_of_S1 = obj.S1_blocknum*(obj.numof_S1_in_one_block);
            
            
            %add distanceweight to weighting
            change_back=sortrows(change_list,2,'descend');
            change_back=change_back(:,[2 1]);
            for i=1:size(change_back,1)
                prelim_list_s1(prelim_list_s1==change_back(i,1)) = change_back(i,2);
            end
            for i=1:max(size(prelim_list_s1))
                obj.timeline.addweight2supersegmentbyid(prelim_list_s1(i),obj.weight.dist)
            end
            obj.timeline.plotsupersegmentweights(obj.signal.bi_stimchan)
            
        end
        
        
        function findS1S2peaks(obj)
            %determins s1 & s2 peaks depending on the weighting
            
            obj.createthreshold_dist;
            
            %try finding outlying blocks - old - is done outside this
            %function in obj.deletesegmentsoutsideblocks 
            detecting_outliers_finished = true;
            exclusion_cnt=0;
            orig_location = diff(obj.S1_blocknum_locs(:,3)); %uses block centers
            reduced_location = diff(obj.S1_blocknum_locs(:,3));
%             while detecting_outliers_finished~=true
%                 num_of_clusters = 2;
%                 if isempty(reduced_location)
%                     cluster_ids = [];
%                 else
%                     [cluster_ids dist_to_centroid] = kmeans(reduced_location,num_of_clusters);
%                 end
%                 outlier_id = cluster_ids~=mode(cluster_ids);
%                 obj.variables.block_deviation_threshold = 2*obj.S1_time; %TODO think about this
%                 
%                 if abs(mean(reduced_location(~outlier_id))-mean(reduced_location(outlier_id)))>=obj.variables.block_deviation_threshold
%                     reduced_location = reduced_location(~outlier_id);
%                     exclusion_cnt=exclusion_cnt+1;
%                 else
%                     fprintf('Finished detecting outliers. %i segment(s) excluded.\n',exclusion_cnt);
%                     detecting_outliers_finished = true;
%                 end
%             end
            
            if exclusion_cnt~=0
                %correct blocknum
                obj.S1_blocknum = obj.S1_blocknum-exclusion_cnt;
                
                %getting from dist back to real vals
                out_ids = [find(~ismember(orig_location,reduced_location)) find(~ismember(orig_location,reduced_location))+1];
                out_vals = orig_location(out_ids);
                %find closest to reduced_loc
                [~,out_id_final] = max(out_vals - mean(reduced_location));
                good_blocknumloc_ids = find([1:max(size(obj.S1_blocknum_locs))~=out_id_final]);
                
                %get s1s in that block
                idlist = obj.gets1ids_inblock_fromblockcenter_id(good_blocknumloc_ids);
                %get all other S1 that are not in idlist
                all_good_ids = obj.timeline.getsupersegids_larger_than_weight(obj.threshold.dist-1);
                badidlist = all_good_ids(~ismember(all_good_ids,idlist(:)));
                
                %get min & max samplerate of those
                %for i=1:size(idlist(:))
                %    minfs(i) = obj.timeline.supersegments{idlist(i)}.min_fs; %- 5 because it makes them special
                %    maxfs(i) = obj.timeline.supersegments{idlist(i)}.max_fs;
                %end
                %min_final=min(minfs)+1; %This cuts of S2 stimuli at the end....
                %max_final=max(minfs)+1;
                
                %delete weighting of bad segments
                for i=1:size(badidlist(:))
                    obj.timeline.supersegments{badidlist(i)}.weight = obj.threshold.pre_dist-5; %- 5 because it makes them special
                end
                
                if obj.settings.ploton
                    obj.timeline.plotsupersegmentweights(obj.signal.bi_stimchan)
                end
            end
            
            
            %delete weights left oft first block
%             outlier = find(~ismember(orig_location,reduced_location));
%             if ~isempty(outlier)
%                 all_ids = 1:max(size(obj.S1_blocknum_locs)); %all ids are rows
%                 used_ids = all_ids(all_ids~=outlier);
%                 minid = find(min(obj.S1_blocknum_locs(used_ids,1))==obj.S1_blocknum_locs(:,1));
%             else %no blocks were deletet, take smallest
%                 minid = find(min(obj.S1_blocknum_locs(:,1))==obj.S1_blocknum_locs(:,1)); %minid is row
%             end
%             %delete weights right of last block
%             cutoff_rightside = max(obj.S1_blocknum_locs(:,2))+ obj.S1_time;
%             for i=1:obj.timeline.numofsupersegments
%                 if obj.timeline.supersegments{i}.max_fs > cutoff_rightside
%                     obj.timeline.supersegments{i}.weight = 0;
%                 end
%             end
            
            %Delete weights outside of blocks
            %TODO why not just count ids backwards from blockcenter?
            %THIS IS DONW WITH A DIFFERENT FUNCTION OUTSIDE OF THIS
            %FUNCTION
%             tmp_ids = obj.gets1ids_inblock_fromblockcenter_id(minid);
%             if ~isempty(tmp_ids)
%                 for i=1:max(size(tmp_ids))
%                     minvals(i) = obj.timeline.supersegments{tmp_ids(i)}.min_fs;
%                 end
%                 startingpoint_fs = floor(min(minvals) - obj.guess_siglength_stim/2); %just as a buffer
%                 
%                 for i=1:obj.timeline.numofsupersegments
%                     if obj.timeline.supersegments{i}.max_fs <= startingpoint_fs
%                         obj.timeline.supersegments{i}.weight = 0;
%                     end
%                 end
%                 
%             end
            
            
            
            %get all good ids that passed the test weight
            all_good_ids = obj.timeline.getsupersegids_larger_than_weight(obj.threshold.dist-1);
            answer='n';
            
            if obj.settings.manualmode==1
                fprintf('Choose threshold.\n')
                f1 = obj.timeline.plotsupersegmentweights(obj.signal.bi_stimchan);
            while ~strcmp(answer,'y')
                figure(f1);
                [~,yval] = ginput(1);
                all_good_ids = obj.timeline.getsupersegids_larger_than_weight(yval-1);
                fprintf('Expected %i points. Found %i points. \n',(obj.S1_blocknum)*(obj.numof_S1_in_one_block+1),numel(all_good_ids))
                
                answer = input('Are you happy with this? [y,n]: ','s');
                
                if isempty(answer) answer='y'; end
            end
            else
                answer = 'y'; 
            end
            num_expected_s1 = obj.numof_S1_in_one_block * obj.S1_blocknum;
            num_expected_s2 = obj.S1_blocknum;
            
            %first try if weighting is sufficiant to detect s1 & s2
            allweights = obj.timeline.getuniqueweights;
            max_weight = max(allweights);
            S1_id_list = obj.timeline.get_segids_by_weight(max_weight);
            
            s2_weight_min = (size(obj.signal.bi_stimchan,2)-1)*(obj.weight.detected_in_chan+obj.weight.hp_detected_in_chan)+(obj.weight.detected_in_stim_chan + obj.weight.hp_detected_in_stim_chan);
            S2_id_list_min = obj.timeline.get_segids_in_weight_window(s2_weight_min-1,max_weight); %-1 to include lower boundary
            if numel(allweights)==1
                s2_weight_max = max_weight;
            else
                s2_weight_max = max(allweights(allweights~=max_weight)); %secondhighest weight
            end
            S2_id_list_max = obj.timeline.get_segids_by_weight(s2_weight_max);
            
            if max(size(S1_id_list)) == num_expected_s1 && max(size(S2_id_list_min)) == num_expected_s2 %if the right amount of peaks for each group was found and weighting doesnt matter
                fprintf('CARE!!! S2 could be wrong, because no check if is really last stimulus of block!\n')
                for i=1:max(size(S1_id_list))
                    obj.timeline.supersegments{S1_id_list(i)}.addsegtype('S1');
                end
                for i=1:max(size(S2_id_list_min))
                    obj.timeline.supersegments{S2_id_list_min(i)}.addsegtype('S2');
                end
                fprintf('Finished!. Found satisfactory S1&S2 peaks through number of peaks.\n')
                obj.timeline.plots1s2(obj.signal.bi_stimchan)
                obj.init.s1s2found = true;
            elseif max(size(S1_id_list)) == num_expected_s1 && max(size(S2_id_list_max)) == num_expected_s2 %if weighting divided peaks perfectly
                for i=1:max(size(S1_id_list))
                    obj.timeline.supersegments{S1_id_list(i)}.addsegtype('S1');
                end
                for i=1:max(size(S2_id_list_max))
                    obj.timeline.supersegments{S2_id_list_max(i)}.addsegtype('S2');
                end
                fprintf('Finished!. Found satisfactory S1&S2 peaks through perfect weighting.\n')
                %for ni=1:9
                ni = obj.stimchan.idpostsplit;
                obj.timeline.plots1s2(obj.signal.bi_stimchan(:,ni))
                %end
                obj.init.s1s2found = true;
            elseif max(size(all_good_ids)) == num_expected_s1 + num_expected_s2 %weighting wrong and
                %reduced_distance = obj.timeline.createreduceddistmat(all_good_ids);
                fprintf('Bad signal or bad algorithms1, trying to reconstruct signal.\n')
                
                %Use counter
                fprintf('Using counter to get S1&S2 - Else use manual method commented in code.\n')
                cnt_s1 = 0;
                %cnt_s2 = 0;
                for i=1:numel(all_good_ids)
                    cnt_s1 = cnt_s1 + 1;
                    if cnt_s1==obj.numof_S1_in_one_block + 1 %meaning s2
                        obj.timeline.supersegments{all_good_ids(i)}.addsegtype('S2');
                        cnt_s1=0;
                    else
                        obj.timeline.supersegments{all_good_ids(i)}.addsegtype('S1');
                    end
                end
                
                %manually select s2
%                 fprintf('Manually select start & stop of S2 stim.\n')
%                 manual_seg = segmentegmtool(obj.signal.bi_split.(obj.stimchan.cath),obj.leads.bi_split.(obj.stimchan.cath));
%                 segments_start = floor(manual_seg(1:2:end,1));
%                 segments_end = floor(manual_seg(2:2:end,1));
%                 
%                 
%                 for i=1:numel(S1_id_list) %set all to S1 and later set some to s2
%                     obj.timeline.supersegments{S1_id_list(i)}.addsegtype('S1');
%                 end
%                 for i=1:numel(all_good_ids)
%                     for j=1:numel(segments_start)
%                         if obj.timeline.supersegments{all_good_ids(i)}.peaks(1)>=segments_start(j) && obj.timeline.supersegments{all_good_ids(i)}.peaks(1)<=segments_end(j)
%                             obj.timeline.supersegments{all_good_ids(i)}.addsegtype('S2');
%                         end
%                     end
%                 end
                
                %obj.reconstructs1s2(obj.numof_S1_in_one_block*obj.S1_time);
                
            else %Try reconstruction of signal
                fprintf('Bad signal or bad algorithms2, trying to reconstruct signal.\n')
                %Create logicline and strfind all real segments and then
                %check if missing, if yes check wich times are big enough
                %to hold segments
                
                obj.reconstructs1s2(obj.numof_S1_in_one_block*obj.S1_time);
                if obj.settings.manualmode==1
                    answer = input('Do you wish to edit this manually? [y,n]: ','s');
                else
                    answer = 'n';
                end
                if isempty(answer) answer='n'; end
                if strcmp(answer,'n')
                    %do nothing for now
                else
                    fprintf('Use spacebar to segment all S1 and S2. Closest supersegment to your manual selection will be used.\n')
                    manual_seg = segmentegmtool(obj.signal.bi_split.(obj.stimchan.cath),obj.leads.bi_split.(obj.stimchan.cath));
                    segments = floor(manual_seg(1:1:end,1));
                    
                    for i=1:obj.timeline.numofsupersegments %resets all supersegments
                        obj.timeline.supersegments{i}.segtype=[];
                    end
                    
                    %go through manual selection and choose closest
                    %supersegment. 
                    for i=1:numel(segments)
                        for j=1:obj.timeline.numofsupersegments
                            diff_tmp(j,i) = abs(obj.timeline.supersegments{j}.peaks(1)-segments(i));
                        end
                    end
                    [~,minidx] = min(diff_tmp,[],1); %saves index of supersegment
                    
                    %Use counter to set s1 and s2 flag
                    cnt_s1 = 0;
                    for i=1:numel(segments)
                        cnt_s1 = cnt_s1 + 1;
                        if cnt_s1==obj.numof_S1_in_one_block + 1 %meaning s2
                            obj.timeline.supersegments{minidx(i)}.addsegtype('S2');
                            cnt_s1=0;
                        else
                            obj.timeline.supersegments{minidx(i)}.addsegtype('S1');
                        end
                    end
                    
                end
                
                
            end
            
            %from knwoledge of blocks & num create perfect sequence
            %expected_sequence_one_block = [zeros(obj.numof_S1_in_one_block,1)' 1 ];
            %expected_sequence_two_blocks = [zeros(obj.numof_S1_in_one_block,1)' 1 zeros(obj.numof_S1_in_one_block,1)'];
            %expected_sequence_full = [zeros(obj.numof_S1_in_one_block,1)' 1 zeros(obj.numof_S1_in_one_block,1)'];
            
            %             for i=1:max(size(S2_id_list)-1)
            %                 abs(obj.timeline.supersegments{S2_id_list(i+1)}.peaks(1)-obj.timeline.supersegments{S2_id_list(i)}.peaks(1))
            %             end
            
            obj.atrial_segments.s1_list = obj.timeline.getS1list_chan_all;
            
        end
        
        function reconstructs1s2(obj,ex)
            %remove this function since it doesnt really work well. Let the
            %programm die in peace
            all_peak_ids_unique = obj.timeline.getsupersegids_larger_than_weight(obj.threshold.pre_dist);
            
            reduced_distance = obj.timeline.createreduceddistmat(all_peak_ids_unique);
            
            mapnew2oldid = [ all_peak_ids_unique;[1:1:max(size(all_peak_ids_unique))] ]'; %col1: supersegid - col2: enum
            change_id = diff(mapnew2oldid,1,2)~=0;
            change_list = mapnew2oldid(change_id,:);
            for i=1:size(change_list,1)
                reduced_distance(reduced_distance(:,1:2)==change_list(i,1)) = change_list(i,2);
            end
            
            %rough reduction
            %find neighbouring ids with dist<2*s1
            neighb_ids = diff([reduced_distance(:,1) reduced_distance(:,2)],1,2)==1;
            neighb_distances = reduced_distance(neighb_ids,:);
            
            %create prelimlist with S1=1
            forward_dist_id = find(neighb_distances(:,3)<=obj.S1_time+obj.threshold.deviation_S1_block_guess & neighb_distances(:,3)>=obj.S1_time-obj.threshold.deviation_S1_block_guess);
            prelim_list_s1 = unique(neighb_distances(forward_dist_id,1:2));
            
            %converting to logical list
            logical_list = zeros(max(size(all_peak_ids_unique)),1);
            logical_list(prelim_list_s1) = 1;
            expected_sequence_one_block = [ones(obj.numof_S1_in_one_block,1)' 0];
            s1_blocks_start = strfind(logical_list',expected_sequence_one_block);
            if max(size(s1_blocks_start)) > 4 %TODO totally random
                obj.init.s1s2found = true;
            end
            
            for i=1:max(size(s1_blocks_start))
                all_s1_ids(i,:) = s1_blocks_start(i):s1_blocks_start(i)+obj.numof_S1_in_one_block;
            end
            all_s2_ids = all_s1_ids(:,end);
            all_s1_ids = all_s1_ids(:,1:end);%end-1
            
            for i=1:max(size(all_s1_ids(:)))
                obj.timeline.supersegments{all_peak_ids_unique(all_s1_ids(i))}.addsegtype('S1');
            end
            for i=1:max(size(all_s2_ids))
                obj.timeline.supersegments{all_peak_ids_unique(all_s2_ids(i))}.addsegtype('S2');
            end
            ids = 1:size(obj.signal.bi_stimchan,2);
            ids = ids(~ismember(ids,[obj.stimchan.idpostsplit,obj.stimchan.idpostsplit+1,obj.stimchan.idpostsplit-1]));
            obj.timeline.plots1s2(obj.signal.bi_stimchan(:,ids))
        end
        
        function [idlist] = gets1ids_inblock_fromblockcenter_id(obj,good_block_ids)
            %rows are good_block_ids
            %columns are supersegids in that block
            
            center_vals = obj.S1_blocknum_locs(good_block_ids);
            
            all_peak_ids_unique = obj.timeline.getsupersegids_larger_than_weight(obj.threshold.pre_dist);
            reduced_distance = obj.timeline.createreduceddistmat(all_peak_ids_unique);
            
            mapnew2oldid = [ all_peak_ids_unique;[1:1:max(size(all_peak_ids_unique))] ]'; %col1: supersegid - col2: enum
            change_id = diff(mapnew2oldid,1,2)~=0;
            change_list = mapnew2oldid(change_id,:);
            for i=1:size(change_list,1)
                reduced_distance(reduced_distance(:,1:2)==change_list(i,1)) = change_list(i,2);
            end
            
            %rough reduction
            %find neighbouring ids with dist<2*s1
            neighb_ids = diff([reduced_distance(:,1) reduced_distance(:,2)],1,2)==1;
            neighb_distances = reduced_distance(neighb_ids,:);
            
            %create prelimlist with S1=1
            forward_dist_id = find(neighb_distances(:,3)<=obj.S1_time+obj.threshold.deviation_S1_block_guess & neighb_distances(:,3)>=obj.S1_time-obj.threshold.deviation_S1_block_guess);
            prelim_list_s1 = unique(neighb_distances(forward_dist_id,1:2));
            
            %converting to logical list
            logical_list = zeros(max(size(all_peak_ids_unique)),1);
            logical_list(prelim_list_s1) = 1;
            expected_sequence_one_block = [ones(obj.numof_S1_in_one_block,1)' 0];
            s1_blocks_start_id = strfind(logical_list',expected_sequence_one_block);
            s1_blocks_end_id = s1_blocks_start_id + obj.numof_S1_in_one_block-1;
            conv2oldid_start = all_peak_ids_unique(s1_blocks_start_id);
            conv2oldid_end = all_peak_ids_unique(s1_blocks_end_id);
            %find s1 blocks closest to center
            for i=1:max(size(s1_blocks_start_id))
                block_boundaries(i,1) = obj.timeline.supersegments{conv2oldid_start(i)}.min_fs;%min
                block_boundaries(i,2) = obj.timeline.supersegments{conv2oldid_end(i)}.max_fs;%min
            end
            %manual method
            %bb = segmentegmtool(obj.signal.bi_stimchan,obj.leads.bi_split.(obj.stimchan.cath));
            
            hitlist = [];
            for i=1:max(size(center_vals))
                for j=1:max(size(s1_blocks_start_id))
                    %find(center_vals(i) >= block_boundaries(:,1) & center_vals(i) <= block_boundaries(:,2))
                    if center_vals(i) >= block_boundaries(j,1) && center_vals(i) <= block_boundaries(j,2)
                        hitlist(i,1) = i;%center_id
                        hitlist(i,2) = block_boundaries(j,1);
                        hitlist(i,3) = block_boundaries(j,2);
                    end
                end
            end
            
            idlist=[];
            for i=1:size(hitlist,1)
                idcnt=0;
                for j=1:size(all_peak_ids_unique,2)
                    if mean(obj.timeline.supersegments{all_peak_ids_unique(j)}.peaks) <= hitlist(i,3) && mean(obj.timeline.supersegments{all_peak_ids_unique(j)}.peaks) >= hitlist(i,2)
                        idcnt=idcnt+1;
                        idlist(i,idcnt) = obj.timeline.supersegments{all_peak_ids_unique(j)}.segment_id;
                    end
                end
            end
            
            
        end
        
        function refine_segments(obj)
            %uses diff instead of nleo for s2 detection
            stimchan_ids = 1:obj.stimchan.size_chan; 
            stimchan_ids(obj.stimchan.exclude_chan) = []; %this option can
            %improve or kill everything
            [~,score,~,~,~,~] = pca(obj.signal.bi_stimchan(:,stimchan_ids));
            sig = score(:,1);
            
            timesegs_s1 = obj.timeline.getsupersegbystring('S1');
            for supsegs=1:numel(timesegs_s1)
                [~,loc]=max(abs(diff(sig(timesegs_s1(supsegs,1).min_fs:timesegs_s1(supsegs,1).max_fs))));
                timesegs_s1(supsegs,1).peaks = loc+timesegs_s1(supsegs,1).min_fs-1;
                x_win_right = obj.threshold.stimuluslength_fs;
                %x_win_right = loc-1; %mirrow left window to right
                timesegs_s1(supsegs,1).max_fs = timesegs_s1(supsegs,1).peaks+x_win_right-1; %create new timewindows
                %plot(timesegs_s1(supsegs,1).peaks,0.1,'*','LineWidth',3)
            end
            obj.timeline.setsupersegbystring('S1',timesegs_s1);
            
            if obj.settings.ploton
                figure
                hold on
                plot(sig) %obj.signal.bi_stimchan(:,1)
            end
            timesegs_s2 = obj.timeline.getsupersegbystring('S2');
            for supsegs=1:numel(timesegs_s2)
                [~,loc]=max(abs(diff(sig(timesegs_s2(supsegs,1).min_fs:timesegs_s2(supsegs,1).max_fs))));
                timesegs_s2(supsegs,1).peaks = loc+timesegs_s2(supsegs,1).min_fs-1;
                x_win_right = obj.threshold.stimuluslength_fs;%(obj.threshold.stimuluslength_fs+obj.threshold.stimuluslength_fs_deviation)/2;
                %x_win_right = loc-1; %mirrow left window to right
                timesegs_s2(supsegs,1).max_fs = timesegs_s2(supsegs,1).peaks+x_win_right-1; %create new timewindows
                if obj.settings.ploton
                    plot(timesegs_s2(supsegs,1).peaks,0.1,'*','LineWidth',3);
                    patch([timesegs_s2(supsegs,1).min_fs,timesegs_s2(supsegs,1).min_fs,timesegs_s2(supsegs,1).max_fs,timesegs_s2(supsegs,1).max_fs],[-ceil(1),ceil(1),ceil(1),-ceil(1)],'red','FaceColor', [255 150 0]/255, 'FaceAlpha', .5, 'EdgeColor', 'none');
                end
                %plot([timesegs_s2(supsegs,1).min_fs:timesegs_s2(supsegs,1).max_fs],sig(timesegs_s2(supsegs,1).min_fs:timesegs_s2(supsegs,1).max_fs),'-k')
                
                %plot([timesegs_s2(supsegs,1).min_fs:timesegs_s2(supsegs,1).max_fs-1],diff(sig(timesegs_s2(supsegs,1).min_fs:timesegs_s2(supsegs,1).max_fs)),'-r')
                
            end
            obj.timeline.setsupersegbystring('S2',timesegs_s2);
            
        end
        
        function gets2times(obj)
            %gets s2 times & refines S1 guess
            obj.timeline.update_numofsupersegments
            s1timelist = 0;
            s2timelist = 0;
            lasts1 = 0;
            s1counter = 0;
            obj.timeline.update_numofsupersegments
            for i=1:obj.timeline.numofsupersegments
                if strcmp(obj.timeline.supersegments{i}.segtype,'S1')
                    [~,max_sig_loc]= max(diff(obj.signal.bi_stimchan((obj.timeline.supersegments{i}.min_fs:obj.timeline.supersegments{i}.max_fs))));
                    max_sig_loc = max_sig_loc + obj.timeline.supersegments{i}.min_fs - 1; %global loc
                    s1timelist(end+1) = max_sig_loc;%mean(obj.timeline.supersegments{i}.peaks);
                    s1counter = s1counter + 1;
                    if s1counter==obj.numof_S1_in_one_block
                        s1counter = 0;  %reset
                        lasts1(end+1) = s1timelist(end);%mean(obj.timeline.supersegments{i}.peaks);
                    end
                elseif strcmp(obj.timeline.supersegments{i}.segtype,'S2')
                    [~,max_sig_loc]= max(diff(obj.signal.bi_stimchan((obj.timeline.supersegments{i}.min_fs:obj.timeline.supersegments{i}.max_fs))));
                    max_sig_loc = max_sig_loc + obj.timeline.supersegments{i}.min_fs - 1;
                    
                    s2timelist(end+1) = max_sig_loc;%mean(obj.timeline.supersegments{i}.peaks);
                end
            end
            s1timelist = sort(s1timelist(2:end),'ascend'); %because of end+1 first value always 0;
            s2timelist = sort(s2timelist(2:end),'ascend');
            lasts1 = lasts1(2:end);
            
            if isempty(s1timelist) || isempty(s2timelist) || isempty(lasts1)
                fprintf('No s1 or s2 found. check dis. \n')
                return
            end
            
            if max(size(lasts1))==max(size(s2timelist))
                timediffs = diff([lasts1;s2timelist],1);
            else
                fprintf('wrong s1s2 detection. found getting s2times\n');
                timediffs=[];
            end
            
            obj.s2_fs = s2timelist;
            obj.s2_steps = timediffs;%fs
            diff_s1 = diff(s1timelist);
            ignore_s1_diff_ids = [0:obj.numof_S1_in_one_block:max(size(diff_s1))];
            ignore_s1_diff_ids = ignore_s1_diff_ids(2:end); %because first val is good and shouldnt be ignored
            diff_s1(ignore_s1_diff_ids) = [];
            
            %refined s1
            %obj.S1_time = mean(diff_s1);
            
            fprintf('S2 times: %.3f \n',obj.s2_steps/obj.samplerate)
            fprintf('S2 timesteps: %.3f\n',abs(diff(obj.s2_steps/obj.samplerate,1,2)))
            %fprintf('refined S1 time: %.3f\n',obj.S1_time)
        end
        
        function getatrialtimewindow_s1(obj)
            
            s1_locs = obj.timeline.getstimuluslocationtable('S1');
            lasts2window = min(obj.s2_steps); %smalles window
            if lasts2window>= 0.9*obj.S1_time %if window > 90% of S1 distance -> to close to next stimulus
                fprintf('Bad S1 window chosen, using half the window\n')
                lasts2window=lasts2window/2;
            end
            
            k=obj.variables.k;%use existing k as upper variable because it was used to detect stims
            s1_val = 2; %we chose END of Stimulation as sttartingpoint for our window
            for i=1:size(obj.signal.bi_stimchan,2)
                cntj = 0;
                for j=1:size(s1_locs,1)
                    cntj=cntj+1;
                    atrial_window_current(i,cntj,:) = [s1_locs(j,s1_val) (s1_locs(j,s1_val)+lasts2window)]; %rows are atrial segments , columns: all samples
                    atrial_window_s1_and_aa(i,cntj,:) = [s1_locs(j,1) (s1_locs(j,1)+lasts2window)]; %rows are atrial segments , columns: all samples
                    
                    signal_segment_aa(i,cntj,:) = obj.signal.bi_stimchan(atrial_window_current(i,cntj,1):atrial_window_current(i,cntj,2),i)'; %chan x s2_id x signal
                    signal_s1_with_aa{i,cntj,:} = obj.signal.bi_stimchan(atrial_window_s1_and_aa(i,cntj,1):atrial_window_s1_and_aa(i,cntj,2),i)';
                end
            end
            
            %             figure
            %             hold on
            %             for i=1:size(signal_segment,1)
            %                 plot(squeeze(signal_segment(i,end-2,:)))
            %             end
            %             hold off
            
            %check if initial NLEO found smth in each sequence else run
            %NLEO %TODO for now just run NLEO
            for i=1:size(obj.signal.bi_stimchan,2) %chan
                cntj = 0;
                for j=1:size(s1_locs,1) %s2
                    cntj=cntj+1;
                    %run NLEO on segment
                    %nleo_filt(i,cntj,:) = nleo( lowpass(squeeze(signal_segment(i,cntj,:)),1,obj.samplerate),obj.samplerate,1,k);
                    nleo_signal_segment(i,cntj,:) = nleo(squeeze(signal_segment_aa(i,cntj,:)),obj.samplerate,1,k);
                    [step(i,cntj,:),activseg(i,cntj,:)] = getActiveSegmentsFromNLEOInShortWindow(squeeze(nleo_signal_segment(i,cntj,:)),obj.samplerate,k);
                    %[step_filt(i,cntj,:),activseg_filt(i,cntj,:)] = getActiveSegmentsFromNLEOInShortWindow(squeeze(nleo_filt(i,cntj,:)),obj.samplerate,k);
                    %activseg{i,cntj} = signal.activseg(i)
                end
            end
            obj.atrial_segments.s1_aa_signal_segment = signal_segment_aa;
            obj.atrial_segments.s1_aa_window_fs = atrial_window_s1_and_aa;
            obj.atrial_segments.s1_aa_window_with_s1_fs = signal_s1_with_aa;
            %obj.atrial_segments.s1_aa_nleo = nleo_signal_segment;
            %obj.atrial_segments.s1_aa_nleo_filt = nleo_filt;
            %obj.atrial_segments.s1_step = step;
            %obj.atrial_segments.s1_activseg = activseg;
            %obj.atrial_segments.s1_step_filt = step_filt;
            %obj.atrial_segments.s1_activseg_filt = activseg_filt;
            
        end
        
        function checkVF(obj)
            s1list = obj.timeline.getS1list_all;
            %s1_list_start_min = repmat(min(s1list,[],1),[obj.stimchan.size_chan 1 1]);
            %s1_list_end_max = repmat(max(s1list,[],1),[obj.stimchan.size_chan 1 1]);
            %s1list(:,:,1) = s1_list_start_min(:,:,1);%replaces starting value with earliest detection for all channels
            %s1list(:,:,2) = s1_list_end_max(:,:,2);
            s1list_start_end = s1list(:,1:2);
            %find s1 overlapping with ventricular far field
            for s1s=1:size(s1list_start_end,1)
                overlap_ventricl_s1 = [];
                win_fs_s1{s1s,1} = [s1list_start_end(s1s,1):s1list_start_end(s1s,2)];
                for ecgs=1:size(obj.ecg.detection,1)
                    ecg_win_fs{ecgs,1} = obj.ecg.detection(ecgs,5):obj.ecg.detection(ecgs,7);
                    winsize_s1=size(win_fs_s1{s1s,1},2);
                    overlap_ventricl_s1(ecgs) = sum(ismember(win_fs_s1{s1s,1}(ceil(winsize_s1/2):end),ecg_win_fs{ecgs,1}))>0; %overlap only if overlap is with righthand half of the s2 window
                end
                obj.ecg.overlap_s1_ventricl(s1s) = sum(overlap_ventricl_s1);
            end
            
            %In order to compare use same beginning time to cut out s2
            %segments.
            s2list = obj.timeline.getS2list_all;
            %s2_list_start_min = repmat(min(s2list,[],1),[obj.stimchan.size_chan 1 1]);
            %s2_list_end_max = repmat(max(s2list,[],1),[obj.stimchan.size_chan 1 1]);
            %s2list(:,:,1) = s2_list_start_min(:,:,1);%replaces starting value with earliest detection for all channels
            %s2list(:,:,2) = s2_list_end_max(:,:,2);
            %s2list_start_end = squeeze(mean(squeeze(s2list(:,:,1:2)),1)); %take mean of all channels
            s2list_start_end = s2list(:,1:2);
            
           
            %find s2 overlapping with ventricular far field
            for s2s=1:size(s2list_start_end,1)
                overlap_ventricl_s2 = [];
                win_fs_s2{s2s,1} = [s2list_start_end(s2s,1):s2list_start_end(s2s,2)];
                for ecgs=1:size(obj.ecg.detection,1)
                    ecg_win_fs{ecgs,1} = obj.ecg.detection(ecgs,5):obj.ecg.detection(ecgs,7);
                    winsize_s2=size(win_fs_s2{s2s,1},2);
                    overlap_ventricl_s2(ecgs) = sum(ismember(win_fs_s2{s2s,1}(ceil(winsize_s2/2):end),ecg_win_fs{ecgs,1}))>0; %overlap only if overlap is with righthand half of the s2 window
                end
                obj.ecg.overlap_s2_ventricl(s2s) = sum(overlap_ventricl_s2);
            end
            
            if obj.settings.ploton
                figure
                hold on
                chan=1;
                plot(obj.ecg.filt.ydata(:,1)) %ecg
                plot(obj.signal.bi_stimchan(:,chan)) %CS
                y_max = max(obj.ecg.filt.ydata(:,chan));
                for i=1:size(obj.ecg.detection,1)
                    patch([obj.ecg.detection(i,5),obj.ecg.detection(i,5),obj.ecg.detection(i,7),obj.ecg.detection(i,7)],[-ceil(y_max),ceil(y_max),ceil(y_max),-ceil(y_max)],'red','FaceColor', [255 150 0]/255, 'FaceAlpha', .5, 'EdgeColor', 'none');
                end
                for ii=1:size(s2list,1)%size(s2list,2)
                    %patch([s2list(chan,ii,1),s2list(chan,ii,1),s2list(chan,ii,2),s2list(chan,ii,2)],[-ceil(y_max),ceil(y_max),ceil(y_max),-ceil(y_max)],'red','FaceColor', [155 50 0]/255, 'FaceAlpha', .5, 'EdgeColor', 'none');
                    patch([s2list(ii,1),s2list(ii,1),s2list(ii,2),s2list(ii,2)],[-ceil(y_max),ceil(y_max),ceil(y_max),-ceil(y_max)],'red','FaceColor', [155 50 0]/255, 'FaceAlpha', .5, 'EdgeColor', 'none');
                end
            end
        end
        
        function getPwavestats(obj)
            s2list = obj.timeline.getS2list_all;
            s2list_start_end = s2list(:,1:2);
            
%             figure
%             hold on
%             chan=1;
%             plot(obj.ecg.filt.ydata(:,1),'r') %ecg
%             plot(obj.signal.bi_stimchan(:,chan),'b') %CS
            
            %Normal detection of p waves doesnt work for stimulation
            %sequence. Use timewindow of Stimulus, check for sequence with
            %nleo and use that as p wave approximation?
            %p = obj.ecg.detection(:,1:3);
            %go through detected P waves
            for pws=1:size(obj.ecg.detection,1)
                p(pws,1) = obj.ecg.detection(pws,1); %ponset
                %check if was already used
            end
%             for s2s=1:size(s2list_start_end,1)
%                 overlap_ventricl_s2 = [];
%                 win_fs_s2{s2s,1} = [s2list_start_end(s2s,1):s2list_start_end(s2s,2)];
%                 for ecgs=1:size(obj.ecg.detection,1)
%                     ecg_win_fs{ecgs,1} = obj.ecg.detection(ecgs,5):obj.ecg.detection(ecgs,7);
%                     winsize_s2=size(win_fs_s2{s2s,1},2);
%                     overlap_ventricl_s2(ecgs) = sum(ismember(win_fs_s2{s2s,1}(ceil(winsize_s2/2):end),ecg_win_fs{ecgs,1}))>0; %overlap only if overlap is with righthand half of the s2 window
%                 end
%                 obj.ecg.overlap_s2_ventricl(s2s) = sum(overlap_ventricl_s2);
%             end
        end
        
        function decideStimTemplate(obj)
            s2list = permute(repmat(obj.timeline.getS2list_all,[1 1 obj.stimchan.size_chan]),[3 1 2]);
            
            %if last window overlaps with VFF ask user if it's ok to use it
            if obj.ecg.overlap_s2_ventricl(end)==1
                fprintf('Last S2 seems to overlap with ventricular farfield.\n')
                alternative_pick = find(obj.ecg.overlap_s2_ventricl==0,1,'last');
                fprintf('Check if ERP has been reached here.\n')
                
                if obj.settings.manualmode==1
                    figure
                    plot(obj.signal.bi_stimchan(:,1))
                    patch([s2list(1,alternative_pick,1),s2list(1,alternative_pick,1),s2list(1,alternative_pick,2),s2list(1,alternative_pick,2)],[-ceil(1),ceil(1),ceil(1),-ceil(1)],'red','FaceColor', [155 50 0]/255, 'FaceAlpha', .5, 'EdgeColor', 'none');
                
                    answer = input('Is the template ok? [y,n]: ','s');
                else
                    answer = 'y';
                end
                %if isempty(answer) answer='y'; end
                if answer=='n'
                    answer = input('Give Block_id of template you want to use [Int]: ','s');
                    s2_block_pick = str2double(answer);
                else
                    s2_block_pick = alternative_pick;
                end
            else
                s2_block_pick = size(s2list,2);
                if obj.settings.manualmode==1
                    answer = input('Using last block as template. Press Enter to continue or input the id you want to use instead [Int]: ','s');
                else
                    answer = '';
                end
                if ~isempty(answer)
                    s2_block_pick = str2double(answer);
                end
            end
            
            obj.variables.s2_block_pick_template = s2_block_pick;
        end
        
        function gets1activity(obj)
            %In order to compare use same beginning time to cut out s1
            %segments.
%             s1list = obj.timeline.getS1list_chan_all;
%             s1_list_start_min = repmat(min(s1list,[],1),[obj.stimchan.size_chan 1 1]);
%             s1_list_end_max = repmat(max(s1list,[],1),[obj.stimchan.size_chan 1 1]);
%             s1list(:,:,1) = s1_list_start_min(:,:,1);%replaces starting value with earliest detection for all channels
%             s1list(:,:,2) = s1_list_end_max(:,:,2);
            s1list = permute(repmat(obj.timeline.getS1list_all,[1 1 obj.stimchan.size_chan]),[3 1 2]);
            
            %In order to compare use same beginning time to cut out s2
            %segments.
%             s2list = obj.timeline.getS2list_chan_all;
%             s2_list_start_min = repmat(min(s2list,[],1),[obj.stimchan.size_chan 1 1]);
%             s2_list_end_max = repmat(max(s2list,[],1),[obj.stimchan.size_chan 1 1]);
%             s2list(:,:,1) = s2_list_start_min(:,:,1);%replaces starting value with earliest detection for all channels
%             s2list(:,:,2) = s2_list_end_max(:,:,2);
            s2list = permute(repmat(obj.timeline.getS2list_all,[1 1 obj.stimchan.size_chan]),[3 1 2]);
            
            lasts2window = min(obj.s2_steps);
            s2list_start_end = squeeze(mean(squeeze(s2list(:,:,1:2)),1)); %take mean of all channels
            
            s2_block_pick = obj.variables.s2_block_pick_template;
            
            for chan=1:obj.stimchan.size_chan
                
                s2window = s2list(chan,s2_block_pick,1):s2list(chan,s2_block_pick,2); %last s2
                processed_sig = obj.signal.bi_stimchan(:,chan);
                
                last_s2_template = squeeze(obj.signal.bi_stimchan(s2window,chan)); %last s2
                
                %think about getting threshold from lonely s2 at end 10% of baseline or
                %somth like that.
                for i=1:size(s1list,2)
                    %fprintf('Iteration: chan: %i s1stim: %i\n',chan,i) %For Inforpurposes
                    idstart(chan,i) = s1list(chan,i,1); %start s1 stim
                    idend(chan,i) = idstart(chan,i) + lasts2window;%obj.atrial_segments.s1_aa_window_fs(chan,i,2); %end of window including atrial activity
                    id_window(chan,i,:) = [idstart(chan,i) idend(chan,i)];
                    stimAAsig{chan,i} = obj.signal.bi_stimchan(idstart(chan,i):idend(chan,i),chan); %window including atrial activity after each S1
                    [~,stimAAsig_as(chan,i)] = getActiveSegmentsFromNLEOInShortWindow(nleo(stimAAsig{chan,i},obj.samplerate,1),obj.samplerate,3);
                    
                    subtracted{chan,i} = filter_stimartifact_matchedfilt(stimAAsig{chan,i},last_s2_template,10,2);
                    processed_sig(idstart(chan,i):idend(chan,i)) = subtracted{chan,i};
                    
                    
                    %kill 80% of stim by setting to baseline
                    percentage_nulled = 100;
                    baseline = mean(obj.signal.bi_stimchan(s1list(chan,end,1)-11:s1list(chan,end,1)-1,chan)); %10 previous samples
                    max_perc_id = ceil(percentage_nulled*numel(stimAAsig{chan,i})/100);
                    ids(chan,i) = idstart(chan,i);
                    if isempty(stimAAsig_as{chan,i}) || numel(stimAAsig_as{chan,i})==1
                        ide(chan,i) = NaN; %This is till end of active segment
                    else
                        ide(chan,i) = idstart(chan,i) + stimAAsig_as{chan,i}(2)-1; %This is till end of active segment
                        processed_sig(ids(chan,i):ide(chan,i)) = baseline;
                    end
                    
                    if ~isnan(ids(chan,i)) && ~isnan(ide(chan,i))
                        smoothwindow = 10; %num of samples ranging into the signal should be smoothed
                        baseline_2smooth = ids(chan,i):ide(chan,i);
                        sig2_smooth = baseline_2smooth(end)+1:baseline_2smooth(end)+smoothwindow;
                        smooth_ids_glob = [baseline_2smooth sig2_smooth];
                        transition{i} = smooth(processed_sig(smooth_ids_glob),20,'loess');
                        
                        processed_sig(smooth_ids_glob) = transition{i};
                    end
                    
                    atrialactivity{chan,i} = processed_sig(idstart(chan,i):idend(chan,i));
                    
                    tmp = nleo(atrialactivity{chan,i},obj.samplerate,1);
                    
                    %                     %retry getting atrial activity segments with processed
                    %                     %signal
                    [~,aa(chan,i)] = getActiveSegmentsFromNLEOInShortWindow(tmp,obj.samplerate,3);
                    %                     for ii=1:numel(aa)
                    %                         if size(aa{ii},1) == 1 && size(aa{ii},2) == 2 %only take properly detected singular atrial activity detect by nleo to calc mean
                    %                             d_time(ii) = diff(aa{ii},1,2);
                    %                         end
                    %                     end
                    %                     mean_aa_samples=round(median(d_time));
                    %
                    %                     for ii=1:size(aa{chan,:},2) %get mean distance for each channel seperately, since distances to source are different
                    %                         if size(aa{ii},1) == 1 && size(aa{ii},2) == 2
                    %                             %s1aa_start(ii) = aa{chan,ii}(end,1);
                    %                             %s1aa_end(ii) = aa{chan,ii}(end,2);
                    %                             s1aa_d(ii) = aa{chan,ii}(end,1);
                    %                         end
                    %                     end
                    %                     mean_dist_s1_aa = round(median(s1aa_d));
                    
                    if size(aa(chan,i),1) >= 1 && size(aa(chan,i),2) == 2
                        s1_list_AA(chan,i,1) = aa{chan,i}(1) + idstart(chan,i) - 1;
                        s1_list_AA(chan,i,2) = aa{chan,i}(2) + idstart(chan,i) - 1;
                    else
                        %s1_list_AA(chan,i,1) = mean_dist_s1_aa + idstart(chan,i) - 1;
                        %s1_list_AA(chan,i,2) = mean_dist_s1_aa + mean_aa_samples + idstart(chan,i) - 1;
                    end
                    [~, nl_loc(chan,i)] = max(tmp);
                    s1_list_AA(chan,i,3) = nl_loc(chan,i) + idstart(chan,i) - 1;
                    
                    
                    %aa_lf{chan,i} = ECG_Low_Filter(atrialactivity{chan,i},obj.samplerate,10);
                    %aa_lf{chan,i} = aa_lf{chan,i}.^2;
                    %gf = fit([1:numel(aa_lf{chan,i})]',aa_lf{chan,i},'gauss1');
                    %nl_loc(chan,i) = gf.b1;
                    %[~, nl_loc(chan,i)] = max(diff(aa_lf{chan,i}));
                    aa_loc_global(chan,i) = nl_loc(chan,i) + idstart(chan,i);
                end
                obj.signal.processed_s1_stimchan(:,chan) = processed_sig;
            end
            
            
            obj.atrial_segments.s1_aa_loc_of_aa = aa_loc_global;
            obj.atrial_segments.s1_aa_sig_without_stim = atrialactivity;
            obj.atrial_segments.s1_aa_window_with_s1_fs = stimAAsig;
            obj.atrial_segments.s1_aa_window_with_s1 = id_window;
            %Plot test
            %
%                         for chan=6
%                             for sigid=1%:obj.numof_S1_in_one_block:size(s1list,2) %try first Atrial Signal of S1 beat in one block
%                                 figure
%                                 hold on
%                                 name = ['chan: ' num2str(chan) ' ' 'stim: ' num2str(sigid)];
%                                 ax1 = subplot(5,1,1);
%                                 plot(stimAAsig{chan,sigid},'Color','b')
%                                 legend('measured')
%                                 ax2 = subplot(5,1,2);
%                                 plot(nleo(stimAAsig{chan,sigid},obj.samplerate,1),'Color','r')
%                                 legend('NLEO measured')
%                                 ax3 = subplot(5,1,3);
%                                 plot(subtracted{chan,sigid},'Color','k')
%                                 legend('subtracted')
%                                 ax4 = subplot(5,1,4);
%                                 plot(atrialactivity{chan,sigid},'Color','g')
%                                 legend('subtracted+nulled')
%                                 ax5 = subplot(5,1,5);
%                                 hold on
%                                 plot(nleo(atrialactivity{chan,sigid},obj.samplerate,1),'Color','r')
%                                 plot(nl_loc(chan,sigid),0,'r-*')
%                                 %check if diff is same as nleo
%                                 [~, diff_loc(chan,sigid)] = max(diff(atrialactivity{chan,sigid}));
%                                 plot(diff_loc(chan,sigid),0,'y-o')
%                                 linkaxes([ax1,ax2,ax3,ax4,ax5],'x');
%                                 legend('NLEO')
%                                 title(name)
%                                 hold off
%                             end
%                         end
            
        end
        
        function detectERP(obj,threshold)
            %detects ERP & missing values
            %and creates a table of good/trustworthy values
            %output:
            %ERP
            %list of trustworthy blocks
            %list of channels to exclude
            
            %get baseline from s2_aa_windows for each atrial activity
            %get amplitude of s1_aa to find general amplitude for the
            %channel & measurement
            
            %get p2p from stim & p2p from s1_aa and compare ratio
            
            s1_perc_aa = obj.signal_properties.s1_aa_p2p./obj.signal_properties.s1_p2p;
            s1_to_aa = median(s1_perc_aa,2); %median per channel
            s2_perc_aa = obj.signal_properties.s2_aa_p2p./obj.signal_properties.s2_p2p;
            
            
            %threshold = obj.signal_properties.s2_window_baseline;
            %threshold = 0.10; %10percent
            
            %second method, just get amplitude of S1 AA and compare
            %s1_to_aa = median(obj.signal_properties.s1_aa_p2p,2);
            %s2_perc_aa = obj.signal_properties.s2_aa_p2p;
            
            %Find values under threshold
            excludetable = zeros(size(s2_perc_aa));
            [excl_r excl_c]= find(s2_perc_aa < threshold*s1_to_aa); %use 10% of ratio of S1 to its aa for S2 to its aa
            
            excludetable(excl_r,excl_c) = 1;
            trusttable = not(excludetable);
            
            numofaa_after_s2_in_all_chan = sum(trusttable);
            
            [km_blocks,km_blocks_idx] = kmeans(numofaa_after_s2_in_all_chan',2,'Start',[min(numofaa_after_s2_in_all_chan);max(numofaa_after_s2_in_all_chan)]); %2 clusters -> stimblocks with signals and without signals
            %count backwards to find first block with signals, since last should have no sig for sure!
            bad_aa_blocks = find(numofaa_after_s2_in_all_chan<=min(km_blocks_idx));
            bad_aa_blocks_logic = logical(zeros(1,size(numofaa_after_s2_in_all_chan,2)));
            bad_aa_blocks_logic(bad_aa_blocks)=1;
            erp_id_guess = find(bad_aa_blocks_logic~=1,1,'last') + 1;%+1 because erp_id is defined as id in wich no stimulus is detected, not as id in wich last stim is detected
            if isempty(erp_id_guess); erp_id_guess=numel(bad_aa_blocks_logic); bad_aa_blocks_logic(:)=0; end
            if erp_id_guess>size(obj.S1_blocknum_locs,1) erp_id_guess=numel(bad_aa_blocks_logic); end
          
            %fprintf('ERP (first block without atrial signals): %d\n',erp_id_guess)
            
            %=========================
            %find channels to exclude, e.g. how many channels next to
            %stimchannel can be ignored
            %find them based on TODO: s1_aa s2_aa and more than one missing
            %=========================
            %[r,~] = find(trusttable(:,bad_aa_blocks_logic)==1);
            %[r,~] = find(excludetable(:,1:erp_id_guess-1)==1);
            %exclude_chan = unique(r);
            %exclude_chan = [];
            
            %obj.stimchan.exclude_chan = exclude_chan;
            s1_aa_nonstim = median(obj.signal_properties.nonstim_s1_aa_p2p,2); %median per channel
            
            %threshold = 0.10; %10percent
            %Find values under threshold
            excludetable_nonstim = zeros(size(obj.signal_properties.nonstim_s2_aa_p2p));
            [excl_r excl_c]= find(obj.signal_properties.nonstim_s2_aa_p2p < threshold*s1_aa_nonstim); %use 10% of ratio of S1 to its aa for S2 to its aa
            
            excludetable_nonstim(excl_r,excl_c) = 1;
            trusttable_nonstim = not(excludetable_nonstim);
            
            numofaa_after_s2_in_all_chan = sum(trusttable);
            
            [km_blocks,km_blocks_idx] = kmeans(numofaa_after_s2_in_all_chan',2,'Start',[min(numofaa_after_s2_in_all_chan);max(numofaa_after_s2_in_all_chan)]); %2 clusters -> stimblocks with signals and without signals
            %count backwards to find first block with signals, since last should have no sig for sure!
            bad_aa_blocks_nonstim = find(numofaa_after_s2_in_all_chan<=min(km_blocks_idx));
            bad_aa_blocks_logic_nonstim = logical(zeros(1,size(numofaa_after_s2_in_all_chan,2)));
            bad_aa_blocks_logic_nonstim(bad_aa_blocks_nonstim)=1;
            erp_id_guess_nonstim = find(bad_aa_blocks_logic_nonstim~=1,1,'last') + 1;%+1 because erp_id is defined as id in wich no stimulus is detected, not as id in wich last stim is detected
            if isempty(erp_id_guess_nonstim); erp_id_guess_nonstim=numel(bad_aa_blocks_logic_nonstim); bad_aa_blocks_logic_nonstim(:)=0; end
            if erp_id_guess_nonstim>size(obj.S1_blocknum_locs,1) erp_id_guess_nonstim=numel(bad_aa_blocks_logic_nonstim); end
            
            
%             %==========OLD WAY, slightly adjusted for paper======
%             %Use non-stimchannel as well since AAs well defined
%             nonstim_sig = ECG_High_Filter(obj.signal.bi_split.(obj.stimchan.non_stimchan),obj.samplerate,5);
%             %nonstim_sig = ECG_Baseline_Removal(obj.signal.bi_split.(obj.stimchan.non_stimchan)(:,3),obj.samplerate,0.01,0.1);
%             
%             s2list = obj.atrial_segments.s2_aa_list; %<--- used clinically
%             
%             s_choice = 1; %clinicaly I used 1, 3 makes more sense
%             s2_win_min = ceil(mean(s2list(:,:,s_choice),1,'omitnan'));
%             s2_delta_cent = ceil(mean(s2list(:,:,3),1,'omitnan')) - s2_win_min;
%             s2_win_max = s2_win_min + obj.s2_steps(end);%max(obj.atrial_segments.s2_aa_window_fs(:,:,2),[],1);
%             
%             %s2list = obj.timeline.getS2list_all;
%             %s2_win_min = ceil(s2list(:,1));
%             %s2_win_max = s2_win_min + obj.s2_steps(end);
%             %TODO
%             %move detection of nonstimchanpeaks to seperate fkt.
%             %use wavelet to detect stimulus
%             
%             %PLotting AA windows to check
%             %             for i=1:size(nonstim_sig,2)
%             %                 figure
%             %                 plot(nonstim_sig(:,i),'Color','b')
%             %                 for j=1:size(s2_win_min,2)
%             %                     hold on
%             %                     patch([s2_win_min(j),s2_win_min(j),s2_win_max(j),s2_win_max(j)],[-ceil(max(abs(nonstim_sig(:,i)))),ceil(max(abs(nonstim_sig(:,i)))),ceil(max(abs(nonstim_sig(:,i)))),-ceil(max(abs(nonstim_sig(:,i))))],'red','FaceColor', [220 20 60]/255, 'FaceAlpha', .5, 'EdgeColor', 'none')
%             %                 end
%             %                 hold off
%             %             end
%             
%             %find good k
%             k_nonstim = obj.variables.k - obj.init.step_k;
%             k_guess_nonstim = false;
%             cnt = 0;
%             remember_best = zeros(1,2);
%             zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0,2,'last'); 
%             while k_guess_nonstim ~= true
%                 cnt = cnt + 1;
%                 k_nonstim = k_nonstim + obj.init.step_k;
%                 
% %                 %run NLEO on the nonstim atrial timewindows
%                  clearvars aa_sig nleo_aa_sig activseg step excludetable_nonstim trusttable_nonstim
%                 for blockid = 1:size(s2_win_min,2)
%                     for chan = 1:size(obj.atrial_segments.nonstim_s2_aa_sig_without_stim,1)
%                         if ~isnan(s2_win_min(blockid))
%                             aa_sig(:,chan,blockid) = obj.atrial_segments.nonstim_s2_aa_sig_without_stim{chan,blockid};
%                             nleo_aa_sig(:,chan,blockid) = nleo(aa_sig(:,chan,blockid),obj.samplerate,1);
%                             [step(chan,blockid,:),~] = getActiveSegmentsFromNLEOInShortWindow(squeeze(nleo_aa_sig(:,chan,blockid)),obj.samplerate,k_nonstim);
%                             as = zci(squeeze(step(chan,blockid,:))-0.95);
%                             if isempty(as)
%                                 activseg{chan,blockid} = [];     
%                             elseif length(as) == 1 || as(2) == size(step,3)
%                                 activseg{chan,blockid} = as(1);
%                             else
%                                 activseg{chan,blockid} = as;
%                             end
%                         end
%                     end
%                 end
%                 
% %                 for blockid = 1:size(s2_win_min,2)
% %                     if ~isnan(s2_win_min(blockid))
% %                         aa_sig(:,:,blockid) = nonstim_sig(s2_win_min(blockid):s2_win_max(blockid),:);
% %                         aa_sig(:,:,blockid) = ECG_Low_Filter(squeeze(aa_sig(:,:,blockid)),obj.samplerate,150);
% %                         nleo_aa_sig(:,:,blockid) = nleo(aa_sig(:,:,blockid),obj.samplerate,1);
% %                         [step(:,blockid,:),activseg(:,blockid,:)] = getActiveSegmentsFromNLEOInShortWindow(squeeze(nleo_aa_sig(:,:,blockid)),obj.samplerate,k_nonstim);
% %                     end
% %                 end
% %                aa_sig = permute(aa_sig,[2 3 1]);%chan x block x signal
%                 excludetable_nonstim = cellfun(@(x) numel(x)<2,activseg); %exclude mean no aa detected
%                 %excludetable_nonstim = cellfun(@(x) numel(x)<1,activseg); %good in theory for sythetic data
%                 %TIPP: for clinical data maybe use amplitude of all found
%                 %signals in nostim and exclude smaller ones
%                 trusttable_nonstim = not(excludetable_nonstim); %detected AA
%                 
%                 if sum(excludetable_nonstim(:)) > remember_best(1,1)
%                     remember_best(1,1) = sum(excludetable_nonstim(:));
%                     remember_best(1,2) = k_nonstim;
%                 end
%                 %if sum(excludetable_nonstim(:)) > numel(excludetable_nonstim)*0.80 %more than 80% of all values are good - debatable if this is a good criteria
%                 if sum(trusttable_nonstim(:)) > numel(excludetable_nonstim)*0.80
%                     k_guess_nonstim = true;
%                 end
%                 if cnt >= obj.threshold.max_itterations_S1_guess*10
%                     break
%                 end
%             end
%             
%             
%             
%             %PLotting AA windows to check
%             %             for i=1:size(nonstim_sig,2)
%             %                 figure
%             %                 plot(nonstim_sig(:,i),'Color','b')
%             %                 for j=1:size(s2_win_min,2)
%             %                     if trusttable_nonstim(i,j)==1
%             %                         hold on
%             %                         patch([activseg{i,j,1}(1)+s2_win_min(j),activseg{i,j,1}(1)+s2_win_min(j),activseg{i,j,1}(2)+s2_win_min(j),activseg{i,j,1}(2)+s2_win_min(j)],[-ceil(max(abs(obj.signal.bi_stimchan(:,i)))),ceil(max(abs(obj.signal.bi_stimchan(:,i)))),ceil(max(abs(obj.signal.bi_stimchan(:,i)))),-ceil(max(abs(obj.signal.bi_stimchan(:,i))))],'red','FaceColor', [220 20 60]/255, 'FaceAlpha', .5, 'EdgeColor', 'none')
%             %                     end
%             %                 end
%             %                 hold off
%             %             end
%             
%             %excludetable_nonstim = cellfun(@isempty,activseg);
%             %trusttable_nonstim = not(excludetable_nonstim);
%             
%             numofaa_after_s2_in_all_chan_nonstim = sum(trusttable_nonstim);
%             
%             %bad_aa_blocks_nonstim = find(numofaa_after_s2_in_all_chan_nonstim<=floor(mode(numofaa_after_s2_in_all_chan_nonstim/2))); %2 el always have stimulus because close to there, and one
%             %bad_aa_blocks_nonstim_logic = logical(zeros(1,size(numofaa_after_s2_in_all_chan_nonstim,2)));
%             %bad_aa_blocks_nonstim_logic(bad_aa_blocks_nonstim)=1;
%             %erp_id_guess_nonstim = find(bad_aa_blocks_nonstim_logic~=1,1,'last') + 1;%+1 because erp_id is defined as id in wich no stimulus is detected, not as id in wich last stim is detected
%             
%             erp_id_guess_nonstim = find(numofaa_after_s2_in_all_chan_nonstim<=1,1,'first');
%             %=================================================
            
            
            
            fprintf('ERP (first block without atrial signals): stimchan: %d nonstimchan: %d \n',erp_id_guess,erp_id_guess_nonstim)
            %[r_nonstim,~] = find(trusttable_nonstim(:,bad_aa_blocks_nonstim_logic)==1);
            
            erp = min([erp_id_guess,erp_id_guess_nonstim]);
            while erp>size(obj.S1_blocknum_locs,1) erp=erp-1; end
            erpblock = obj.S1_blocknum_locs(erp,1:2);
            
            if obj.settings.ploton
                figure
                hold on
                plot(obj.signal.bi_split.(obj.stimchan.non_stimchan)(:,1))
                y_max = max(obj.signal.bi_split.(obj.stimchan.non_stimchan)(:,1));
                patch([erpblock(1,1),erpblock(1,1),erpblock(1,2),erpblock(1,2)],[-ceil(y_max),ceil(y_max),ceil(y_max),-ceil(y_max)],'red','FaceColor', [255 150 0]/255, 'FaceAlpha', .5, 'EdgeColor', 'none')
                title('ERP in NON-stimulationcatheter')
                figure
                hold on
                allchanstim = [1:obj.stimchan.size_chan];
                allchanstim_shift = circshift(allchanstim,obj.stimchan.size_chan-obj.stimchan.idpostsplit);
                goodchan = allchanstim_shift(end-floor(numel(allchanstim_shift)/2));
                %plot(obj.signal.bi_split.(obj.stimchan.cath)(:,goodchan))
                plot(ECG_Baseline_Removal(obj.signal.bi_split.(obj.stimchan.cath)(:,goodchan),obj.samplerate,0.01,0.1))
                patch([erpblock(1,1),erpblock(1,1),erpblock(1,2),erpblock(1,2)],[-ceil(y_max),ceil(y_max),ceil(y_max),-ceil(y_max)],'red','FaceColor', [255 150 0]/255, 'FaceAlpha', .5, 'EdgeColor', 'none')
                title('ERP in stimulationcatheter')
            end
            
            if obj.settings.manualmode==1
                answer = input('Do you wish to edit this manually? [y,n]: ','s');
            else
                answer = 'n';
            end
            if isempty(answer) answer='n'; end
            %Here there is ERP detection from nonstim channels
            if answer=='n'
                if erp_id_guess_nonstim==erp_id_guess | isempty(erp_id_guess_nonstim)
                    obj.erp_id = erp_id_guess;
                    obj.bad_aa_blocks = bad_aa_blocks_logic;
                else
                    %obj.erp_id = erp_id_guess;
                    obj.erp_id = erp_id_guess_nonstim;
                    obj.bad_aa_blocks = bad_aa_blocks_logic;
                end
            else
                obj.erp_id = input('Input your ERP. (It is defined as the first Block WITHOUT an atrial answer to the stimulus: ');
                answer = input('Input ids of bad blocks seperated by , : ','s');
                if ~isempty(answer)
                    conv = regexp(answer,',','split');
                    obj.bad_aa_blocks = cellfun(@str2double, conv);
                else
                    obj.bad_aa_blocks = [];
                end
                
            end
            %Here a manual correction method in case it doesnt work
            %sometimes add this somehow
            
            %             %ask user for ERP
            %             sigplot.title=sprintf('EGM bip');
            %             sigplot.sig{1}.leads=obj.leads.in;
            %             sigplot.sig{1}.signals=obj.signal.in_bi;
            %             sigplot.fs=obj.samplerate;
            %             sigplot.sig{1}.color='b';
            %             plotMultiChannelEGM(sigplot);
            %
            %             %inp = input('What is ERP? Input id of first Block without atrial signal: ');
            %             %obj.erp_id = inp;
            
            
            
        end
        
        function getatrialtimewindow_s2(obj)
            %gets timewindow(s) based on the previously determined
            %stimulus signal (s2):
            %
            %-> s2_aa_signal_segment:
            %signalsegment starts from end of stimulus (s2_end) and ends after
            %the atrial activity ends
            %
            %better would be to introduce window from start of s2 till end
            %of window, but for now this happens in "obj.getatrialactivity"
            %TODO
            
            s2_locs = obj.timeline.getstimuluslocationtable('S2');
            %s2_timesegs = obj.timeline.getsupersegbystring('S2');
            lasts2window = min(obj.s2_steps); %in samples
            fprintf('Rethink this thouroughly!!\n')
            %maybe a  little less than a second
            %keyboard
            k=obj.variables.k; %use existing k as upper variable because it was used to detect stims
            
            s2_val = 2; %2 because 1=start and 2=end of segment
            for i=1:size(obj.signal.bi_stimchan,2)
                cntj = 0;
                for j=1:size(s2_locs,1)
                    cntj=cntj+1;
                    atrial_window_current(i,cntj,:) = [s2_locs(j,s2_val) (s2_locs(j,s2_val)+lasts2window)]; %rows are atrial segments , columns: all samples
                    signal_segment(i,cntj,:) = obj.signal.bi_stimchan(atrial_window_current(i,cntj,1):atrial_window_current(i,cntj,2),i)'; %chan x s2_id x signal
                end
            end
            
            %             figure
            %             hold on
            %             for i=1:size(signal_segment,1)
            %                 plot(squeeze(signal_segment(i,end-2,:)))
            %             end
            %             hold off
            
            %check if initial NLEO found smth in each sequence else run
            %NLEO %TODO for now just run NLEO
            for i=1:size(obj.signal.bi_stimchan,2) %chan
                cntj = 0;
                for j=1:size(s2_locs,1) %s2
                    cntj=cntj+1;
                    %run NLEO on segment
                    %nleo_filt(i,cntj,:) = nleo( lowpass(squeeze(signal_segment(i,cntj,:)),1,obj.samplerate),obj.samplerate,1,k);
                    nleo_signal_segment(i,cntj,:) = nleo(squeeze(signal_segment(i,cntj,:)),obj.samplerate,1,k);
                    [~,activseg(i,cntj,:)] = getActiveSegmentsFromNLEOInShortWindow(squeeze(nleo_signal_segment(i,cntj,:)),obj.samplerate,k);
                    %[step_filt(i,cntj,:),activseg_filt(i,cntj,:)] = getActiveSegmentsFromNLEOInShortWindow(squeeze(nleo_filt(i,cntj,:)),obj.samplerate,k);
                    %activseg{i,cntj} = signal.activseg(i)
                end
            end
            %obj.atrial_segments.s2_aa_signal_segment = signal_segment;
            %obj.atrial_segments.s2_aa_window_fs = atrial_window_current;
            %obj.atrial_segments.s2_aa_nleo = nleo_signal_segment;
            %obj.atrial_segments.s2_aa_nleo_filt = nleo_filt;
            %obj.atrial_segments.s2_activseg = activseg;
            %obj.atrial_segments.s2_step_filt = step_filt;
            %obj.atrial_segments.s2_activseg_filt = activseg_filt;
        end
        
        function getatrialactivity(obj)
            %Evaluation of previously obtained atrial activity window
            
            %Take window with S2 Stimulus
            %remove stimulus
            %save new AA segment & save filtered signal
            
            s2list_old = permute(repmat(obj.timeline.getS2list_all,[ 1 1 10]),[3 1 2]);
            s2list = obj.timeline.getS2list_chan_all;
            %origin_s2 = squeeze(s2list(obj.stimchan.idpostsplit,:,:));%test
            
            %In order to compare use same beginning time to cut out s2
            %segments.
            s2_list_start_min = repmat(min(s2list,[],1),[obj.stimchan.size_chan 1 1]);
            s2_list_end_max = repmat(max(s2list,[],1),[obj.stimchan.size_chan 1 1]); %TODO CHANGED THE FOLLLOWING
            %s2list(:,:,1) = s2_list_start_min(:,:,1);%replaces starting value with earliest detection for all channels
            %s2list(:,:,2) = s2_list_end_max(:,:,2);
            
            lasts2window = min(obj.s2_steps);
            %use minimum NLEO end of all channels for the first few
            %pacingtrains since atrial activity close
            %             replace = 1:3; %replace first till third found NLEO with min NLEO
            %             r = min(s2list(:,replace,2),[],1);
            %             for i=1:size(s2list,1)
            %                 for j=1:size(replace)
            %                     s2list(i,replace(j),2) = r(j);
            %                 end
            %             end
            s2_block_pick = obj.variables.s2_block_pick_template;
            
            for chan=1:obj.stimchan.size_chan
                
                %s2window = s2list(end,1):s2list(end,2);
                s2window = s2list(chan,s2_block_pick,1):s2list(chan,s2_block_pick,2); %last s2
                
                %if numel(s2window)==1 %in case s2list_all has zero entries
                %    s2window = s2list_old(end,1):s2list_old(end,2);
                %end
                processed_sig(:,chan) = obj.signal.bi_stimchan(:,chan);
                
                %maxnum_of_stim = 2;
                %think about getting threshold from lonely s2 at end 10% of baseline or
                %somth like that.
                for i=1:size(s2list,2)
                    
                    last_s2_template = squeeze(obj.signal.bi_stimchan(s2window,chan)); %last s2
                    %get only stimulus with nleo
                    %[~,actseg]=getActiveSegmentsFromNLEOInShortWindow(nleo(last_s2_template,obj.samplerate,1),obj.samplerate,.5);
                    %last_s2_template(1:actseg{1,1}(1))=0;
                    %last_s2_template(actseg{1,1}(2):end)=0;
                    
                    idstart(chan,i) = s2list(chan,i,1); %start s2 stim
                    if idstart(chan,i)==0 %in case nothing was found
                        idstart(chan,i) = s2list_old(i,1); 
                    end
                    
                    idend(chan,i) = s2list(chan,i,2)+lasts2window;%obj.atrial_segments.s2_aa_window_fs(chan,i,2); %end of window including atrial activity
                    id_window(chan,i,:) = [idstart(chan,i) idend(chan,i)];
                    %fprintf('Chan: %i,Sig: %i\n',chan,i)%incase loop fails
                    stimAAsig{chan,i} = obj.signal.bi_stimchan(idstart(chan,i):idend(chan,i),chan); %s2 windows after each block
%                     k=3;
%                     [~,stimAAsig_as(chan,i)] = getActiveSegmentsFromNLEOInShortWindow(nleo(stimAAsig{chan,i},obj.samplerate,1),obj.samplerate,k);
%                     while numel(stimAAsig_as{chan,i})<2 
%                         [~,stimAAsig_as(chan,i)] = getActiveSegmentsFromNLEOInShortWindow(nleo(stimAAsig{chan,i},obj.samplerate,1),obj.samplerate,k);
%                         k = k-obj.init.step_k;
%                         if k <= 0
%                             %stimAAsig_as{chan,i} = [1 stimAAsig_as{chan,i}]; %Quickfix
%                             break;
%                         end
%                     end
                    stimAAsig_as{chan,i} = squeeze(s2list(chan,i,1:2));
                    subtracted{chan,i} = filter_stimartifact_matchedfilt(stimAAsig{chan,i},last_s2_template,10,1);
                    
                    processed_sig(idstart(chan,i):idend(chan,i),chan) = subtracted{chan,i};
                    
                    
%                   figure
%                   hold on
%                   plot(stimAAsig{chan,i})
%                   %plot(zeropad_stim)
%                   plot(subtracted{chan,i})
                     
%                    kill 80% of stim by setting to baseline
%                     percentage_nulled = 10;
%                     baseline = mean(obj.signal.bi_stimchan(s2list(chan,end,1)-11:s2list(chan,end,1)-1,chan)); %10 previous samples
%                     max_perc_id = ceil(percentage_nulled*numel(stimAAsig{chan,i})/100);
%                     processed_sig(s2list(chan,i,1):s2list(chan,i,1)+max_perc_id-1) = baseline;
%                     peakcoice_diff = 1;
%                     pks_sort_diff(peakcoice_diff,2) = max_perc_id;
%                     %processed_sig(s2list(chan,i,1):s2list(chan,i,2)) = baseline;
%                     ids(chan,i) = idstart(chan,i);
%                     ide(chan,i) = ids(chan,i) + max_perc_id;

                     ids(chan,i) = idstart(chan,i);
                     if isempty(stimAAsig_as{chan,i}) || numel(stimAAsig_as{chan,i})==1
                         ide(chan,i) = NaN; %This is till end of active segment
                         sig_wav = method_wavelet(stimAAsig{chan,i}','WaveletShape','bior1.5','samplerate',obj.samplerate,'ThresholdFactor',1,'MinimumInaktivityLength',10,'MinimumAktivityLength',15,'Postprocessing','On');
                         [~,stimAAsig_as(chan,i)] = getActiveSegmentsFromNLEOInShortWindow(sig_wav',obj.samplerate,obj.variables.k);%0.1 -> obj.variables.k
                         if isempty(stimAAsig_as{chan,i}) || numel(stimAAsig_as{chan,i})==1
                             ide(chan,i) = NaN;
                         else
                             diff_sig = abs(diff(stimAAsig{chan,i})); %worked well without abs()
                             %windowstim = 0.025*obj.samplerate;%obj.signal_properties.mean_signallength_fs;
                             %diff_sig = abs(diff(stimAAsig{chan,i}(1:1+windowstim))); %worked well without abs()
                             [s_diff_val, s_diff_loc] = findpeaks( diff_sig );
                             pks_sort_diff = sortrows([s_diff_val s_diff_loc],'descend');
                             buffer=5;%buffer=5;
                             peakcoice_diff = 1; %1 with buffer is better than 2 without buffer!
                             ide(chan,i) = idstart(chan,i) + pks_sort_diff(peakcoice_diff,2)+buffer; %ide si from start till stimulus
                             baseline = linspace(processed_sig(ids(chan,i)-1,chan),processed_sig(ide(chan,i),chan), ide(chan,i)-ids(chan,i)+2);
                             processed_sig(ids(chan,i)-1:ide(chan,i),chan) = baseline;
                         end
                     else
                         %ide(chan,i) = idstart(chan,i) + stimAAsig_as{chan,i}(2)-1; %This is till end of active segment
                         %processed_sig(ids(chan,i):ide(chan,i)) = baseline;
                        
                        diff_sig = abs(diff(stimAAsig{chan,i})); %worked well without abs()
                        %windowstim = 0.025*obj.samplerate;%obj.signal_properties.mean_signallength_fs;
                        %diff_sig = abs(diff(stimAAsig{chan,i}(1:1+windowstim))); %%doesnt work well!!
                        %[~,diff_loc]=max(diff_sig);
                        [s_diff_val, s_diff_loc] = findpeaks( diff_sig );
                        pks_sort_diff = sortrows([s_diff_val s_diff_loc],'descend');
                        buffer=5;%buffer=5;
                        peakcoice_diff = 1; %1 with buffer is better than 2 without buffer!
                        ide(chan,i) = idstart(chan,i) + pks_sort_diff(peakcoice_diff,2)+buffer; %ide si from start till stimulus
                        baseline = linspace(processed_sig(ids(chan,i)-1,chan),processed_sig(ide(chan,i),chan), ide(chan,i)-ids(chan,i)+2);
                        processed_sig(ids(chan,i)-1:ide(chan,i),chan) = baseline;
                    end
%                     
%                     smoothwindow = 80; %num of samples ranging into the signal should be smoothed
%                     smoothpoints = 10; %navX 20
%                     baseline_2smooth = ids(chan,i):ide(chan,i);
%                     sig2_smooth = baseline_2smooth(end)+1:baseline_2smooth(end)+smoothwindow;
%                     smooth_ids_glob = [baseline_2smooth sig2_smooth];
%                     if sum(isnan(smooth_ids_glob))>=1
%                         transition{i} = [];
%                     else
%                         transition{i} = smooth(processed_sig(smooth_ids_glob,chan),smoothpoints,'loess');
%                     end
%                     
%                     processed_sig(smooth_ids_glob,chan) = transition{i};
                    
                    atrialactivity{chan,i} = processed_sig(idstart(chan,i):idend(chan,i),chan);
                    
                    tmp = nleo(atrialactivity{chan,i},obj.samplerate,1);
                    
                    %retry getting atrial activity segments with processed
                    %signal
                    seg_std = 0.2 ;%2;%0.2 was originally
                    [~,aa(chan,i)] = getActiveSegmentsFromNLEOInShortWindow(tmp,obj.samplerate,seg_std);
                    fprintf('The problem is here. Make this more sophisticated.\n')
                    %[~,aa(chan,i)] = getActiveSegmentsfromNLEOdiff(tmp);
                    if numel(aa{chan,i})>=2 %TODO changed this
                        s2_list_AA(chan,i,1) = aa{chan,i}(1) + idstart(chan,i) - 2;
                        s2_list_AA(chan,i,2) = aa{chan,i}(2) + idstart(chan,i) - 2;
                        [~, nl_loc(chan,i)] = max(tmp);
                        s2_list_AA(chan,i,3) = nl_loc(chan,i) + idstart(chan,i) - 2;
                        %aa_lf{chan,i} = ECG_Low_Filter(atrialactivity{chan,i},obj.samplerate,10);
                        %aa_lf{chan,i} = aa_lf{chan,i}.^2;
                        %gf = fit([1:numel(aa_lf{chan,i})]',aa_lf{chan,i},'gauss1');
                        %nl_loc(chan,i) = gf.b1;
                        %[~, nl_loc(chan,i)] = max(diff(aa_lf{chan,i}));
                        aa_loc_global(chan,i) = nl_loc(chan,i) + idstart(chan,i) - 2;
                        %aa_s2_dist(chan,i) = abs(aa_loc_global(chan,i) - origin_s2(i,3));
                        aa_s2_dist(chan,i) = nl_loc(chan,i)-pks_sort_diff(peakcoice_diff,2);
                    else
                        s2_list_AA(chan,i,1) = nan;
                        s2_list_AA(chan,i,2) = nan;
                        s2_list_AA(chan,i,3) = nan;
                        aa_loc_global(chan,i) = nan;
                        aa_s2_dist(chan,i) = nan;
                    end
                    
                    
                end
            end
            
            if obj.settings.ploton
                sigplot.title=sprintf('S2 stim + atrial activity');
                sigplot.sig{1}.leads= [repmat({'processed'},numel(obj.leads.bi_split.(obj.stimchan.cath)),1) ; obj.leads.bi_split.(obj.stimchan.cath)];%obj.data.cathdatauni{1,1}.leads;
                sigplot.sig{1}.signals= [processed_sig obj.signal.bi_stimchan];%obj.data.cathdatauni{1,1}.ydata;
                sigplot.fs=obj.samplerate;
                sigplot.sig{1}.color='b';
                plotMultiChannelEGM(sigplot);
            end
            
            answer=[];
            %answer=input('Do you want to manually correct this? (y/n): ','s');
            if strcmp(answer,'y')
                %TODO manual segmentation
                fprintf('TODO\n')
            else
                %save
                obj.signal.processed_s2_stimchan = processed_sig;
                obj.atrial_segments.s2_aa_sig_without_stim = atrialactivity;
                obj.atrial_segments.s2_aa_list = s2_list_AA; %atrial activity found without stimulussig
                obj.atrial_segments.s2_aa_loc_of_aa = aa_loc_global;
                obj.atrial_segments.s2_aa_dist_fs = aa_s2_dist;
                
                obj.atrial_segments.s2_aa_window_with_s2_fs = stimAAsig;
                obj.atrial_segments.s2_aa_window_with_s2 = id_window;
                
                %PLotting S2 windows to check
%                 for i=1:size(s2_list_AA,1)
%                     name = ['chan: ' num2str(obj.stimchan.el1) ' ' 'stim: ' num2str(i)];
%                     figure
%                     plot(obj.signal.bi_stimchan(:,i),'Color','b')
%                     for j=1:size(s2_list_AA,2)
%                         hold on
%                         patch([s2_list_AA(i,j,1),s2_list_AA(i,j,1),s2_list_AA(i,j,2),s2_list_AA(i,j,2)],[-ceil(max(abs(obj.signal.bi_stimchan(:,i)))),ceil(max(abs(obj.signal.bi_stimchan(:,i)))),ceil(max(abs(obj.signal.bi_stimchan(:,i)))),-ceil(max(abs(obj.signal.bi_stimchan(:,i))))],'red','FaceColor', [220 20 60]/255, 'FaceAlpha', .5, 'EdgeColor', 'none')
%                     end
%                     hold off
%                     title(name)
%                 end
            end
            
            %DO NOT DELETE THIS
            %Checks how well eleimination of stimulussignal worked
            
%             for chan = 3%1:obj.stimchan.size_chan
%                 for sigid = 1%:3%size(s2list,2)%3
%                     f = figure;
%                     hold on
%                     name = ['chan: ' num2str(chan) ' ' 'stim: ' num2str(sigid)];
%                     sgtitle(name)
%                     ax1 = subplot(4,1,1);
%                     plot(stimAAsig{chan,sigid},'Color','b','LineWidth',1.5)
%                     set(gca,'Ytick',[]);
%                     set(gca,'Xtick',[]);
%                     %legend('raw signal')
%                     %ax2 = subplot(5,1,2);
%                     %plot(nleo(stimAAsig{chan,sigid},obj.samplerate,1),'Color','r','LineWidth',1)
%                     %set(gca,'Yticklabel',[]);
%                     %legend('NLEO raw signal')
%                     ax3 = subplot(4,1,2);
%                     plot(subtracted{chan,sigid},'Color','k','LineWidth',1.5)
%                     set(gca,'Ytick',[]);
%                     set(gca,'Xtick',[]);
%                     new_xlim = [min(subtracted{chan,sigid}) max(subtracted{chan,sigid})];
%                     %legend('matched filter (MF)')
%                     ax4 = subplot(4,1,3);
%                     plot(atrialactivity{chan,sigid},'Color','k','LineWidth',1.5)
%                     xlim(new_xlim);
%                     set(gca,'Ytick',[]);
%                     set(gca,'Xtick',[]);
%                     %legend('MF + linear interpolated baseline')
%                     ax5 = subplot(4,1,4);
%                     hold on
%                     nlstim = nleo(stimAAsig{chan,sigid},obj.samplerate,1);
%                     plot(nlstim,'-.r','LineWidth',1)
%                     nlaa = nleo(atrialactivity{chan,sigid},obj.samplerate,1);
%                     nlaa = max(nlstim)/max(nlaa) * nlaa;
%                     plot(nlaa,'-.r','LineWidth',1)
%                     plot(nl_loc(chan,sigid),0,'r-*','LineWidth',1)
%                     set(gca,'Yticklabel',[]);
%                     set(gca,'Ytick',[]);
%                     %legend('NLEO processed signal')
%                     %check if diff is same as nleo
%                     [~, diff_loc(chan,sigid)] = max(diff(atrialactivity{chan,sigid}));
%                     %plot(diff_loc(chan,sigid),0,'y-o')
%                     linkaxes([ax1,ax3,ax4,ax5],'x');
%                     xlim([0 200]);
%                     box on
% %                   save_name = 'stimdel';
% %                   savefig(f,['Pics/' save_name '.fig'])
% %                   print(f,['Pics/' save_name],'-dpng','-r400')
%                 end
%            end
%             
            
            %Plotting for EMBC Poster
            %used 4-2 bip and 4-1 cal_bip
            %             for chan = 4%3:obj.stimchan.size_chan
            %                 for sigid = 2
            %
            %                     %figure
            %                     %set(gcf,'Position',[1000 1000 850 260])
            %                     %hold on
            %                     %name = ['chan: ' num2str(chan) ' ' 'stim: ' num2str(sigid)];
            %                     %sgtitle(name)
            %
            %                     %ax1 = subplot(3,1,1);
            %                     figure
            %                     set(gcf,'Position',[1000 1000 850 260])
            %                     hold on
            %                     plot(stimAAsig{chan,sigid},'Color',[60/255 80/255 126/255],'LineWidth',3)
            %                     plot(nleo(stimAAsig{chan,sigid},obj.samplerate,1)*max(stimAAsig{chan,sigid})/max(nleo(stimAAsig{chan,sigid},obj.samplerate,1)),'--r','LineWidth',3)
            %                     box on
            %                     set(gca,'linewidth',2)
            %                     set(gca,'FontSize',25)
            %                     xlim([0 300])
            %                     ylim([-18 12])
            %                     hold off
            %                     %legend('original signal')
            %                     %ax2 = subplot(5,1,2);
            %                     %plot(nleo(stimAAsig{chan,sigid},obj.samplerate,1),'Color','r','LineWidth',3)
            %                     %legend('NLEO')
            %                     %ax3 = subplot(3,1,2);
            %                     figure
            %                     set(gcf,'Position',[1000 1000 850 260])
            %                     plot(subtracted{chan,sigid},'-k','LineWidth',3)
            %                     set(gca,'linewidth',2)
            %                     set(gca,'FontSize',25)
            %                     box on
            %                     xlim([0 300])
            %                     ylim([-18 12])
            %                     %legend('subtracted')
            %                     %ax4 = subplot(3,1,3);
            %                     figure
            %                     set(gcf,'Position',[1000 1000 850 260])
            %                     hold on
            %                     plot(atrialactivity{chan,sigid},'Color',[90/255 120/255 70/255],'LineWidth',3)
            %                     plot(2*nleo(atrialactivity{chan,sigid},obj.samplerate,1)*max(atrialactivity{chan,sigid})/max(nleo(atrialactivity{chan,sigid},obj.samplerate,1)),'--r','LineWidth',3)
            %                     set(gca,'linewidth',2)
            %                     set(gca,'FontSize',25)
            %                     box on
            %                     xlim([0 300])
            %                     ylim([-18 12])
            %                     hold off
            %                     %legend('subtracted+nulled')
            %                     %ax5 = subplot(5,1,5);
            %                     %hold on
            %                     %plot(nleo(atrialactivity{chan,sigid},obj.samplerate,1),'Color','r','LineWidth',3)
            %                     %legend('NLEO')
            %                     %linkaxes([ax1,ax2,ax3,ax4,ax5],'x');
            %                     %linkaxes([ax1,ax3,ax4],'x');
            %                 end
            %             end
            
            %============================
            %testing wavelet
            %============================
            %             x = 3;%ceil(log2(obj.samplerate/2/50));
            %
            %             for chan=1:obj.stimchan.size_chan
            %                 for i=1:size(atrialactivity,2)
            %                     na{chan,i} = atrialactivity{chan,i};
            %                     n = numel(atrialactivity{chan,i});
            %                     while mod(n/2^x,1)~=0
            %                         na{chan,i} = [0; na{chan,i}];
            %                         n = numel(na{chan,i});
            %                     end
            %                     [~,swd{chan,i}]=swt(na{chan,i},x,'bior1.5');
            %                     swd2{chan,i} = abs(swd{chan,i}).^2;
            %                     swd2lp{chan,i} = ECG_Low_Filter(swd2{chan,i},obj.samplerate,1);
            %
            %                     gf = fit([1:numel(swd2lp{chan,i}(end,:))]',swd2lp{chan,i}(end,:)','gauss1');
            %                     atrialactivity_wavelet{chan,i} = gf.b1;
            %                 end
            %             end
            %
            %
            %             for chan=1:9
            %                 figure
            %                 name= [num2str(obj.stimchan.el1) '-' num2str(chan)];
            %                     title(name)
            %                 hold on
            %                 for i=1:8
            %                     dd(chan,i) = (atrialactivity_wavelet{chan,i}+idstart(chan,i)-1-origin_s2(i,2));
            %                 end
            %                 plot(mean(stimchan_data(:,chan))*obj.samplerate./dd(chan,:))
            %                 hold off
            %             end
            %
            %
            %
            %             lvl = 2;
            %             f = obj.samplerate/(2*x);
            %             for chan=1
            %                 for i=1:9
            %                     figure
            %                     ax1 = subplot(4,1,1);
            %                     plot(obj.signal.bi_stimchan(idstart(chan,i):idend(chan,i)))
            %                     ax2 = subplot(4,1,2);
            %                     plot(atrialactivity{chan,i})
            %                     ax3 = subplot(4,1,3);
            %                     plot(swd{chan,i}(lvl,:))
            %                     ax4 = subplot(4,1,4);
            %                     hold on
            %                     plot(swd2lp{chan,i}(lvl,:))
            %                     gf = fit([1:numel(swd2lp{chan,i}(end,:))]',swd2lp{chan,i}(end,:)','gauss1');
            %                     plot(gf)
            %                     hold off
            %                     %findpeaks(swd{chan,i}(lvl,:),obj.samplerate/lvl,'MinPeakHeight',0.1*max(swd{chan,i}(lvl,:)))
            %                     %plot()
            %                     %plot(nleo(swd{chan,i}(lvl,:),obj.samplerate,0))
            %                     name = [num2str(chan) '-' num2str(i)];
            %                     title(name)
            %                     linkaxes([ax1,ax2,ax3,ax4],'x');
            %                 end
            %             end
            
            
            
            
            
            
            
            %check 3 steps
            %             figure
            %             hold on
            %             plot(stimAAsig{chan,i},'Color','b')
            %             plot(zeropad_stim,'Color','k')
            %             plot(subtracted{chan,i},'Color','g')
            %             plot(atrialactivity{chan,i},'Color','r')
            %             hold off
            %             legend('measured','s2template','subtracted','subtracted+nulled')
            %
            %Check whole channel
            %             figure
            %             ax1 = subplot(2,1,1);
            %             plot(obj.signal.bi_stimchan(:,chan))
            %             ax2 = subplot(2,1,2);
            %             plot(processed_sig)
            %             linkaxes([ax1,ax2],'xy');
            
            
            
%             i=21; %signal
%             j=8; %channel
%             
%             figure
%             ax1 = subplot(3,1,1);
%             hold on
%             plot(stimAAsig{j,i},'Color','b')
%             a = getActiveSegmentsFromNLEOInShortWindow(nleo(stimAAsig{j,i},obj.samplerate,1),obj.samplerate,3);
%             patch([find(a==1,1,'first'),find(a==1,1,'first'),find(a==1,1,'last'),find(a==1,1,'last')],[-ceil(max(abs(stimAAsig{j,i}))),ceil(max(abs(stimAAsig{j,i}))),ceil(max(abs(stimAAsig{j,i}))),-ceil(max(abs(stimAAsig{j,i})))],'red','FaceColor', [220 20 60]/255, 'FaceAlpha', .5, 'EdgeColor', 'none')
%             
%             legend('measured','active segment')
%             hold off
%             ax2 = subplot(3,1,2);
%             plot(subtracted{j,i},'Color','k')
%             legend('subtracted')
%             ax3 = subplot(3,1,3);
%             plot(atrialactivity{j,i},'Color','g')
%             legend('subtracted+baselined')
            
            %ax4 = subplot(4,1,4);
            %hold on
            %plot(stimAAsig{j,i},'Color','b')
            %a = getActiveSegmentsFromNLEOInShortWindow(nleo(atrialactivity{j,i},obj.samplerate,1),obj.samplerate,3);
            %patch([find(a==1,1,'first'),find(a==1,1,'first'),find(a==1,1,'last'),find(a==1,1,'last')],[-ceil(max(abs(atrialactivity{j,i}))),ceil(max(abs(atrialactivity{j,i}))),ceil(max(abs(atrialactivity{j,i}))),-ceil(max(abs(atrialactivity{j,i})))],'red','FaceColor', [220 20 60]/255, 'FaceAlpha', .5, 'EdgeColor', 'none')
            %hold off
            %linkaxes([ax1,ax2,ax3,ax4],'x');
            %linkaxes([ax1,ax2,ax3],'x');
            %legend('measured','active segment')
            
            %
            %             for i=1:size(aa_s2_dist,1)
            %                figure
            %                plot(20./aa_s2_dist(i,:)*obj.samplerate)
            %             end
            
            
            
            
            
            
            
            %Matched Filter approach (shorter code):
            % b = flipud(stimsig);
            % y = filter(b,1,stimAAsig{j,i});
            % n = 1:length(y);
            % match = find(y==max(y));
            % figure
            % plot(n,y,'b', n(match), y(match), 'ro');
            % %thresh = 1.1;
            % %norm = b.'*b;
            % %matches = n(y>thresh*norm);
            % figure
            % plot(n,y,'b', n(matches), y(matches), 'ro');
            
            
            
            
            
        end
        
        function combine_segmented_sig(obj)
            %combines signal with removed s1 with the signal with removed
            %s2 stimuli.
            %Uses signal with removed s1 as base and replaces s2 with
            %removed s2 segment
            sig_tmp = obj.signal.processed_s1_stimchan;
            for chan=1:obj.stimchan.size_chan
                for s2blocks = 1:size(obj.timeline.getS2list_chan_all,2)
                    idstart = obj.atrial_segments.s2_aa_window_with_s2(chan,s2blocks,1);
                    idend = obj.atrial_segments.s2_aa_window_with_s2(chan,s2blocks,2);
                    sig_tmp(idstart:idend,chan) = obj.signal.processed_s2_stimchan(idstart:idend,chan);
                end
            end
            obj.signal.processed_all_stimchan = sig_tmp;
        end
        
        function getatrialactivity_unipolar(obj)
            %use atrial activity window from s2list found in bipolar signal
            %to detect steepest slope in unipolar signal
            
            %s2list_stimchan = obj.atrial_segments.s2_aa_list;
            %s2list_nonstimchan = obj.atrial_segments.nonstim_s2_list_AA;
            s2_stim_list = obj.timeline.getS2list_all;
            
            if strcmpi(obj.stimchan.cath,'cs')
                cs_list = obj.atrial_segments.s2_aa_list; %atrial activity after s2 found after filtering out stimulus
                spi_list = obj.atrial_segments.nonstim_s2_list_AA;
            else
                cs_list = obj.atrial_segments.nonstim_s2_list_AA;
                spi_list = obj.atrial_segments.s2_aa_list;
            end
            
            if ~isempty(obj.signal.uni_split.CS)
                for chan=1:size(obj.signal.uni_split.CS,2)
                    for s2block = 1:size(cs_list,2)
                        %retry finding stimulus
                        strt_stim = round(s2_stim_list(s2block,1));
                        endt_stim = round(s2_stim_list(s2block,2));
                        if isnan(strt_stim) || isnan(endt_stim)
                            s2_uni_fs_cs(chan,s2block) = nan;
                        else
                            [~,maxt] = max(obj.signal.uni_split_prefilt.CS(strt_stim:endt_stim,chan));
                            [~,mint] = min(obj.signal.uni_split_prefilt.CS(strt_stim:endt_stim,chan));
                            med = median(obj.signal.uni_split_prefilt.CS(strt_stim:endt_stim,chan));
                            trailing_end = round(0.006*obj.samplerate); %depending on stimgenerator 4-6ms trailing end
                            if maxt>mint %first deflection down then up
                                depol_id=find(diff(obj.signal.uni_split_prefilt.CS(strt_stim:mint+strt_stim-1,chan))>med,1,'last');
                                repol_id=find(diff(obj.signal.uni_split_prefilt.CS(maxt+strt_stim-1:maxt+strt_stim-1+trailing_end,chan))<med,1,'first') + maxt-1;
                            else
                                repol_id=find(diff(obj.signal.uni_split_prefilt.CS(strt_stim:mint+strt_stim-1,chan))>med,1,'last');
                                depol_id=find(diff(obj.signal.uni_split_prefilt.CS(maxt+strt_stim-1:maxt+strt_stim-1+trailing_end,chan))<med,1,'first') + maxt-1;
                            end
                        end
                        if isempty(repol_id) && ~isempty(depol_id)%chose manual window TODO find other algorithm for this
                            repol_id = depol_id + trailing_end;
                            [~,mtdid] = max(diff(obj.signal.uni_split_prefilt.CS(depol_id+strt_stim-1:repol_id+strt_stim-1,chan)));
                        elseif isempty(depol_id) && ~isempty(repol_id)
                            depol_id = repol_id - trailing_end;
                            [~,mtdid] = max(diff(obj.signal.uni_split_prefilt.CS(depol_id+strt_stim-1:repol_id+strt_stim-1,chan)));
                        elseif isempty(depol_id) && isempty(repol_id)
                            depol_id = nan;
                            repol_id = nan;
                            mtdid = nan;
                        else
                            [~,mtdid] = max(diff(obj.signal.uni_split_prefilt.CS(depol_id+strt_stim-1:repol_id+strt_stim-1,chan)));
                            if isempty(mtdid) mtdid=nan; end
                        end
                        %fprintf('Chan %i, Block: %i\n',chan, s2block)
                        %TODO: CONTINUE HERE:
                        %figure
%hold on
%plot(obj.signal.uni_split_prefilt.CS(strt_stim:endt_stim,chan))
%plot(movmean(obj.signal.uni_split_prefilt.CS(strt_stim:endt_stim,chan).^2,10))
                        s2_uni_fs_cs(chan,s2block,:) = [depol_id repol_id mtdid + depol_id] + strt_stim - 1;
                        %detect steepest  of atrial answer
                        strt = cs_list(round(chan/2),s2block,1);
                        endt = cs_list(round(chan/2),s2block,2);
                        if isnan(strt) || isnan(endt)
                            s2peak_uni_fs_cs(chan,s2block) = nan;
                        else
                            [~,mtdid] = max(diff(obj.signal.uni_split_prefilt.CS(strt:endt,chan)));
                            s2peak_uni_fs_cs(chan,s2block) = mtdid + strt;
                        end
                    end
                end
            end
            if ~isempty(obj.signal.uni_split.SPI)
                for chan=1:size(obj.signal.uni_split.SPI,2)
                    for s2block = 1:size(spi_list,2)
                        %retry finding stimulus
                        strt_stim = round(s2_stim_list(s2block,1));
                        endt_stim = round(s2_stim_list(s2block,2));
                        if isnan(strt_stim) || isnan(endt_stim)
                            s2_uni_fs_spi(chan,s2block) = nan;
                        else
                            [~,maxt] = max(obj.signal.uni_split_prefilt.SPI(strt_stim:endt_stim,chan));
                            [~,mint] = min(obj.signal.uni_split_prefilt.SPI(strt_stim:endt_stim,chan));
                            med = median(obj.signal.uni_split_prefilt.SPI(strt_stim:endt_stim,chan));
                            trailing_end = round(0.006*obj.samplerate); %depending on stimgenerator 4-6ms trailing end
                            if maxt>mint %first deflection down then up
                                depol_id=find(diff(obj.signal.uni_split_prefilt.SPI(strt_stim:mint+strt_stim-1,chan))>med,1,'last');
                                repol_id=find(diff(obj.signal.uni_split_prefilt.SPI(maxt+strt_stim-1:maxt+strt_stim-1+trailing_end,chan))<med,1,'first') + maxt-1;
                            else
                                depol_id=find(diff(obj.signal.uni_split_prefilt.SPI(strt_stim:mint+strt_stim-1,chan))>med,1,'last');
                                repol_id=find(diff(obj.signal.uni_split_prefilt.SPI(maxt+strt_stim-1:maxt+strt_stim-1+trailing_end,chan))<med,1,'first') + maxt-1;
                            end
                        end
                        if isempty(repol_id) %chose manual window TODO find other algorithm for this
                            repol_id = depol_id + trailing_end;
                        end
                        [~,mtdid] = max(diff(obj.signal.uni_split_prefilt.SPI(depol_id+strt_stim-1:repol_id+strt_stim-1,chan)));
                        if isempty(mtdid)
                            mtdid = nan;
                        end
                        if isempty(depol_id)
                            depol_id = nan;
                        end
                        uni_acitvity_id = [depol_id repol_id mtdid+depol_id] + strt_stim - 1;
                        if ~isempty(uni_acitvity_id) && size(uni_acitvity_id,2)==3 
                            s2_uni_fs_spi(chan,s2block,:) = uni_acitvity_id;
                        else
                            s2_uni_fs_spi(chan,s2block,:) = [nan nan nan];
                        end
                        %detect steepest slope
                        strt = spi_list(round(chan/2),s2block,1);
                        endt = spi_list(round(chan/2),s2block,2);
                        if isnan(strt) || isnan(endt)
                            s2peak_uni_fs_spi(chan,s2block) = nan;
                        else
                            [~,mtdid] = max(diff(obj.signal.uni_split_prefilt.SPI(strt:endt,chan)));
                            s2peak_uni_fs_spi(chan,s2block) = mtdid + strt;
                        end
                        
                    end
                end
            end
            
            if strcmpi(obj.stimchan.cath,'cs')
                obj.atrial_segments.s2_list_unipolar = s2_uni_fs_cs; %stimulusdetection in unipolar chan
                obj.atrial_segments.s2_aa_unipolar_peak = s2peak_uni_fs_cs;
                obj.atrial_segments.nonstim_s2_unipolar_list = s2_uni_fs_spi;
                obj.atrial_segments.nonstim_s2_aa_unipolar_peak = s2peak_uni_fs_spi;
            else
                obj.atrial_segments.s2_list_unipolar = s2_uni_fs_spi;
                obj.atrial_segments.s2_aa_unipolar_peak = s2peak_uni_fs_spi;
                obj.atrial_segments.nonstim_s2_unipolar_list = s2_uni_fs_cs;
                obj.atrial_segments.nonstim_s2_aa_unipolar_peak = s2peak_uni_fs_cs;
            end
%             %diff(obj.signal.uni_split_prefilt.SPI);
%             for chan=1:size(obj.atrial_segments.s2_aa_list,1)
%                 for s2block = 1:size(obj.atrial_segments.s2_aa_list,2)
%                     strt = s2list_stimchan(chan,s2block,1);
%                     endt = s2list_stimchan(chan,s2block,2);
%                     if isnan(strt) || isnan(endt)
%                         s2peak_uni_fs_cs(chan,s2block) = nan;
%                     else
%                         [~,mtdid] = max(diff(obj.signal.uni_split_prefilt.SPI(strt:endt,chan)));
%                         s2peak_uni_fs_cs(chan,s2block) = mtdid + strt;
%                     end
%                 end
%             end
            
            
            
        end
        
        function getatrialamplitudes(obj)
            %Signal has been segmentetd completly and now we use processed
            %signals and detect the amplitudes of the atrial activities
            %from them
            
            numofsamples_for_baseline_estimation = 20;
            %Amps for atrial activity after S1 Stimuli
            for chan=1:obj.stimchan.size_chan
                for sigid = 1:size(obj.atrial_segments.s1_aa_sig_without_stim,2)
                    baseline_estimate_s1_window(chan,sigid) = mean(obj.atrial_segments.s1_aa_sig_without_stim{chan,sigid}(end-numofsamples_for_baseline_estimation:end)); %take last 20 samples
                    %baseline_estimate_s1_window(chan,sigid) = 0;
                    
                    s1_amp_pos(chan,sigid) = abs(max(obj.atrial_segments.s1_aa_window_with_s1_fs{chan,sigid}) - baseline_estimate_s1_window(chan,sigid)); %TODO
                    s1_amp_neg(chan,sigid) = abs(baseline_estimate_s1_window(chan,sigid) - min(obj.atrial_segments.s1_aa_window_with_s1_fs{chan,sigid}));
                    s1_p2p(chan,sigid) = abs(max(obj.atrial_segments.s1_aa_window_with_s1_fs{chan,sigid}) - min(obj.atrial_segments.s1_aa_window_with_s1_fs{chan,sigid}));
                    
                    s1_aa_amp_pos(chan,sigid) = abs(max(obj.atrial_segments.s1_aa_sig_without_stim{chan,sigid}) - baseline_estimate_s1_window(chan,sigid)); %TODO
                    s1_aa_amp_neg(chan,sigid) = abs(baseline_estimate_s1_window(chan,sigid) - min(obj.atrial_segments.s1_aa_sig_without_stim{chan,sigid}));
                    s1_aa_p2p(chan,sigid) = abs(max(obj.atrial_segments.s1_aa_sig_without_stim{chan,sigid}) - min(obj.atrial_segments.s1_aa_sig_without_stim{chan,sigid}));
                end
            end
            
            
            %Amps for atrial activity after S2 Stimuli
            for chan=1:obj.stimchan.size_chan
                for sigid = 1:size(obj.atrial_segments.s2_aa_sig_without_stim,2)
                    baseline_estimate_s2_window(chan,sigid) = mean(obj.atrial_segments.s2_aa_sig_without_stim{chan,sigid}(end-numofsamples_for_baseline_estimation:end)); %take last 20 samples
                    %baseline_estimate_s2_window(chan,sigid) = 0; %maybe makes more sense to hp filter so that baseline is gone
                    
                    s2_amp_pos(chan,sigid) = abs(max(obj.atrial_segments.s2_aa_window_with_s2_fs{chan,sigid}) - baseline_estimate_s2_window(chan,sigid)); %TODO
                    s2_amp_neg(chan,sigid) = abs(baseline_estimate_s2_window(chan,sigid) - min(obj.atrial_segments.s2_aa_window_with_s2_fs{chan,sigid}));
                    s2_p2p(chan,sigid) = abs(max(obj.atrial_segments.s2_aa_window_with_s2_fs{chan,sigid}) - min(obj.atrial_segments.s2_aa_window_with_s2_fs{chan,sigid}));
                    
                    s2_aa_amp_pos(chan,sigid) = abs(max(obj.atrial_segments.s2_aa_sig_without_stim{chan,sigid}) - baseline_estimate_s2_window(chan,sigid)); %TODO
                    s2_aa_amp_neg(chan,sigid) = abs(baseline_estimate_s2_window(chan,sigid) - min(obj.atrial_segments.s2_aa_sig_without_stim{chan,sigid}));
                    s2_aa_p2p(chan,sigid) = abs(max(obj.atrial_segments.s2_aa_sig_without_stim{chan,sigid}) - min(obj.atrial_segments.s2_aa_sig_without_stim{chan,sigid}));
                end
            end
            
            
            obj.signal_properties.s1_window_baseline = baseline_estimate_s1_window;
            obj.signal_properties.s2_window_baseline = baseline_estimate_s2_window;
            
            obj.signal_properties.s1_p2p = s1_p2p;
            obj.signal_properties.s1_aa_p2p = s1_aa_p2p;
            obj.signal_properties.s1_aa_amp_pos = s1_aa_amp_pos;
            obj.signal_properties.s1_aa_amp_neg = s1_aa_amp_neg;
            obj.signal_properties.s2_p2p = s2_p2p;
            obj.signal_properties.s2_aa_p2p = s2_aa_p2p;
            obj.signal_properties.s2_aa_amp_pos = s2_aa_amp_pos;
            obj.signal_properties.s2_aa_amp_neg = s2_aa_amp_neg;
            
            %Plottest
            %             figure
            %             for chan=1:obj.stimchan.size_chan
            %                 hold on
            %                 plot(flipud(fliplr(s2_aa_p2p(chan,:))),'LineWidth',3)
            %                 leg{chan} = ['Chan: ' num2str(chan)];
            %             end
            %             legend([leg])
            %             hold off
            
        end
        
        function gets1activity_nonstimchan(obj)
            
            s1list = permute(repmat(obj.timeline.getS1list_all,[1 1 size(obj.signal.bi_split.(obj.stimchan.non_stimchan),2)]),[3 1 2]);
            s2list = permute(repmat(obj.timeline.getS2list_all,[1 1 size(obj.signal.bi_split.(obj.stimchan.non_stimchan),2)]),[3 1 2]);
            
            lasts2window = min(obj.s2_steps);
            s2_block_pick = obj.variables.s2_block_pick_template;
            
            for chan=1:size(obj.signal.bi_split.(obj.stimchan.non_stimchan),2)
                
                s2window = s2list(chan,s2_block_pick,1):s2list(chan,s2_block_pick,2); %last s2
                processed_sig = obj.signal.bi_split.(obj.stimchan.non_stimchan)(:,chan);
                
                last_s2_template = squeeze(obj.signal.bi_split.(obj.stimchan.non_stimchan)(s2window,chan)); %last s2
                
                %think about getting threshold from lonely s2 at end 10% of baseline or
                %somth like that.
                for i=1:size(s1list,2)
                    %fprintf('Iteration: chan: %i s1stim: %i\n',chan,i) %For Inforpurposes
                    idstart(chan,i) = s1list(chan,i,1); %start s1 stim
                    idend(chan,i) = idstart(chan,i) + lasts2window;%obj.atrial_segments.s1_aa_window_fs(chan,i,2); %end of window including atrial activity
                    id_window(chan,i,:) = [idstart(chan,i) idend(chan,i)];
                    stimAAsig{chan,i} = obj.signal.bi_split.(obj.stimchan.non_stimchan)(idstart(chan,i):idend(chan,i),chan); %window including atrial activity after each S1
                    [~,stimAAsig_as(chan,i)] = getActiveSegmentsFromNLEOInShortWindow(nleo(stimAAsig{chan,i},obj.samplerate,1),obj.samplerate,3);
                    
                    subtracted{chan,i} = filter_stimartifact_matchedfilt(stimAAsig{chan,i},last_s2_template,10,2);
                    processed_sig(idstart(chan,i):idend(chan,i)) = subtracted{chan,i};
                    
                    
                    %kill 80% of stim by setting to baseline
                    percentage_nulled = 80;
                    baseline = mean(obj.signal.bi_split.(obj.stimchan.non_stimchan)(s1list(chan,end,1)-11:s1list(chan,end,1)-1,chan)); %10 previous samples
                    max_perc_id = ceil(percentage_nulled*numel(stimAAsig{chan,i})/100);
                    ids(chan,i) = idstart(chan,i);
                    if isempty(stimAAsig_as{chan,i}) || numel(stimAAsig_as{chan,i})==1
                        ide(chan,i) = NaN; %This is till end of active segment
                    else
                        ide(chan,i) = idstart(chan,i) + stimAAsig_as{chan,i}(2)-1; %This is till end of active segment
                        processed_sig(ids(chan,i):ide(chan,i)) = baseline;
                    end
                    
                    if ~isnan(ids(chan,i)) && ~isnan(ide(chan,i))
                        smoothwindow = 10; %num of samples ranging into the signal should be smoothed
                        baseline_2smooth = ids(chan,i):ide(chan,i);
                        sig2_smooth = baseline_2smooth(end)+1:baseline_2smooth(end)+smoothwindow;
                        smooth_ids_glob = [baseline_2smooth sig2_smooth];
                        transition{i} = smooth(processed_sig(smooth_ids_glob),20,'loess');
                        
                        processed_sig(smooth_ids_glob) = transition{i};
                    end
                    
                    atrialactivity{chan,i} = processed_sig(idstart(chan,i):idend(chan,i));
                    tmp = nleo(atrialactivity{chan,i},obj.samplerate,1);
                    [~,aa(chan,i)] = getActiveSegmentsFromNLEOInShortWindow(tmp,obj.samplerate,3);
                    
                    if size(aa(chan,i),1) >= 1 && size(aa(chan,i),2) == 2
                        s1_list_AA(chan,i,1) = aa{chan,i}(1) + idstart(chan,i) - 1;
                        s1_list_AA(chan,i,2) = aa{chan,i}(2) + idstart(chan,i) - 1;
                    else
                        %s1_list_AA(chan,i,1) = mean_dist_s1_aa + idstart(chan,i) - 1;
                        %s1_list_AA(chan,i,2) = mean_dist_s1_aa + mean_aa_samples + idstart(chan,i) - 1;
                    end
                    [~, nl_loc(chan,i)] = max(tmp);
                    s1_list_AA(chan,i,3) = nl_loc(chan,i) + idstart(chan,i) - 1;
                    aa_loc_global(chan,i) = nl_loc(chan,i) + idstart(chan,i);
                end
                obj.signal.processed_s1_nonstimchan(:,chan) = processed_sig;
            end
            
            
            obj.atrial_segments.nonstim_s1_aa_loc_of_aa = aa_loc_global;
            obj.atrial_segments.nonstim_s1_aa_sig_without_stim = atrialactivity;
            obj.atrial_segments.nonstim_s1_aa_window_with_s1_fs = stimAAsig;
            obj.atrial_segments.nonstim_s1_aa_window_with_s1 = id_window;
            
%             s1list = permute(repmat(obj.timeline.getS1list_all,[1 1 obj.stimchan.size_chan]),[3 1 2]);
%             
%             tmp = nleo(atrialactivity{chan,i},obj.samplerate,1);
%             [~,aa(chan,i)] = getActiveSegmentsFromNLEOInShortWindow(tmp,obj.samplerate,3);
%             
%             
%             for chan=1:obj.stimchan.size_chan
%                 for i=1:size(s1list,2)
%                     obj.atrial_segments.s1_aa_window_with_s1
%                 end
%             end
%             
%             obj.atrial_segments.nonstim_s1_aa_loc_of_aa %global time of peak center
%             obj.atrial_segments.nonstim_s1_list_AA
        end
        
        function getatrialactivity_nonstimchan(obj)
            
            %getnon-stim channel
            nonstim_sig = obj.signal.bi_split.(obj.stimchan.non_stimchan);
            
            %get atrial window
            
            s2_locs = obj.timeline.getstimuluslocationtable('S2');
            %origin_s2 = squeeze(s2_locs(obj.stimchan.idpostsplit,:));
            %s2list = obj.atrial_segments.s2_aa_list; %Somehow this value is
            %fail.....
            timewindow = round(0.5*obj.samplerate); %Heartbeat is ~0.8s -> P wave shorter
            s2_win_min = min(s2_locs(:,3),[],2);%min(obj.atrial_segments.s2_aa_window_fs(:,:,1),[],1);
            s2_win_max = s2_win_min + timewindow; %this is beginning of next block! %max(obj.atrial_segments.s2_aa_window_fs(:,:,2),[],1);
            fprintf('think about this \n')
            %keyboard
            %             figure
            %             for chan=1:size(nonstim_sig,2)
            %                 for block=1:obj.S1_blocknum
            %                     hold on
            %                     plot(nonstim_sig(s2_win_min(block):s2_win_max(block),chan));
            %                     hold off
            %                 end
            %             end
            
            for chan = 1:size(nonstim_sig,2)
                for blockid = 1:max(size(s2_win_min))
                    aa_sig{chan,blockid} = nonstim_sig(s2_win_min(blockid):s2_win_max(blockid),chan);
                    nleo_aa_sig{chan,blockid} = nleo(aa_sig{chan,blockid},obj.samplerate,1);
                    [~,activseg(chan,blockid)] = getActiveSegmentsFromNLEOInShortWindow(nleo_aa_sig{chan,blockid},obj.samplerate,3);
                    if numel(activseg{chan,blockid})>=2
                        s2_list_AA(chan,blockid,1) = activseg{chan,blockid}(1) + s2_win_min(blockid) - 1;
                        s2_list_AA(chan,blockid,2) = activseg{chan,blockid}(2) + s2_win_min(blockid) - 1;
                        s2_list_AA(chan,blockid,2) = activseg{chan,blockid}(2) + s2_win_min(blockid) - 1;
                        [~, nl_loc(chan,blockid)] = max(nleo_aa_sig{chan,blockid});
                        aa_loc_global(chan,blockid) = nl_loc(chan,blockid) + s2_win_min(blockid) - 1; %start s2 stim
                        s2_list_AA(chan,blockid,3) = aa_loc_global(chan,blockid);
                    else
                        s2_list_AA(chan,blockid,1) = nan;
                        s2_list_AA(chan,blockid,2) = nan;
                        s2_list_AA(chan,blockid,3) = nan;
                        aa_loc_global(chan,blockid) = nan;
                        aa_s2_dist(chan,blockid) = nan;
                    end
                    
                    aa_s2_dist(chan,blockid) = abs(aa_loc_global(chan,blockid) - s2_locs(blockid,3));
                end
            end
            
            
            
            %             figure
            %             for chan=1:size(nonstim_sig,2)
            %                 for blockid=5%1:obj.S1_blocknum
            %                     figure
            %                     hold on
            %                     plot(aa_sig{chan,blockid});
            %                     hold off
            %                 end
            %             end
            
            
            
            %             for chan = 1:4%3:obj.stimchan.size_chan
            %                 for sigid = 1:3
            %                     figure
            %                     hold on
            %                     name = ['chan: ' num2str(chan) ' ' 'stim: ' num2str(sigid)];
            %                     sgtitle(name)
            %                     ax1 = subplot(2,1,1);
            %                     plot(aa_sig{chan,sigid},'Color','b')
            %                     legend('measured')
            %                     ax2 = subplot(2,1,2);
            %                     hold on
            %                     plot(nleo(aa_sig{chan,sigid},obj.samplerate,1),'Color','r')
            %                     plot(nl_loc(chan,sigid),max(nleo(aa_sig{chan,sigid},obj.samplerate,1)),'x')
            %                     hold off
            %                     legend('NLEO measured')
            %                     linkaxes([ax1,ax2],'x');
            %                     legend('NLEO')
            %                 end
            %             end
            
            excludetable_nonstim = cellfun(@isempty,activseg);
            trusttable_nonstim = not(excludetable_nonstim);
            
            numofaa_after_s2_in_all_chan_nonstim = sum(trusttable_nonstim);
            
            obj.atrial_segments.nonstim_s2_aa_sig_without_stim = aa_sig;
            obj.atrial_segments.nonstim_s2_aa_loc_of_aa = aa_loc_global;
            obj.atrial_segments.nonstim_s2_aa_dist_fs = aa_s2_dist;
            obj.atrial_segments.nonstim_s2_list_AA = s2_list_AA; %atrial activity found without stimulussig
            
            %obj.atrial_segments.nonstim_s2_aa_window_with_s2_fs = stimAAsig;
            %obj.atrial_segments.nonstim_s2_aa_window_with_s2 = id_window;
            
        end
        
        function getatrialamplitudes_nonstimchan(obj)
            %getnon-stim channel
            nonstim_sig = obj.signal.bi_split.(obj.stimchan.non_stimchan);
            
            %get atrial window
            
            s2_locs = obj.timeline.getstimuluslocationtable('S2');
            %origin_s2 = squeeze(s2_locs(obj.stimchan.idpostsplit,:,:));
            %s2list = obj.atrial_segments.s2_aa_list; %Somehow this value is
            %fail.....
            s2_win_min = min(s2_locs(:,3),[],2);%min(obj.atrial_segments.s2_aa_window_fs(:,:,1),[],1);
            s2_win_max = s2_win_min + obj.s2_steps(end); %this is beginning of next block! %max(obj.atrial_segments.s2_aa_window_fs(:,:,2),[],1);
            
            %             figure
            %             for chan=1:size(nonstim_sig,2)
            %                 for block=1:obj.S1_blocknum
            %                     hold on
            %                     plot(nonstim_sig(s2_win_min(block):s2_win_max(block),chan));
            %                     hold off
            %                 end
            %             end
            
            
            %Amps for atrial activity after S1 Stimuli
            numofsamples_for_baseline_estimation = 20;
            for chan=1:size(nonstim_sig,2)
                for sigid=1:size(obj.atrial_segments.s1_list,2) %try first Atrial Signal of S1 beat in one block
                    baseline_estimate_s2_window(chan,sigid) = mean(obj.atrial_segments.nonstim_s1_aa_sig_without_stim{chan,sigid}(end-numofsamples_for_baseline_estimation:end)); %take last 20 samples
                    s1_aa_p2p_nonstim(chan,sigid) = max(obj.atrial_segments.nonstim_s1_aa_window_with_s1_fs{chan,sigid}) - min(obj.atrial_segments.nonstim_s1_aa_window_with_s1_fs{chan,sigid});
                end
            end
            obj.signal_properties.nonstim_s1_aa_p2p = s1_aa_p2p_nonstim;
            
            %Amps for atrial activity after S2 Stimuli
            numofsamples_for_baseline_estimation = 20;
            for chan=1:size(nonstim_sig,2)
                for sigid=1:size(obj.atrial_segments.s2_aa_list,2) %try first Atrial Signal of S1 beat in one block
                    baseline_estimate_s2_window(chan,sigid) = mean(obj.atrial_segments.nonstim_s2_aa_sig_without_stim{chan,sigid}(end-numofsamples_for_baseline_estimation:end)); %take last 20 samples
                    s2_aa_p2p(chan,sigid) = max(obj.atrial_segments.nonstim_s2_aa_sig_without_stim{chan,sigid}) - min(obj.atrial_segments.nonstim_s2_aa_sig_without_stim{chan,sigid});
                end
            end
            obj.signal_properties.nonstim_s2_aa_p2p = s2_aa_p2p;
        end
        
        function getlocations(obj)
            if ~isempty(obj.map) && ~isempty(obj.locations)
                obj.locations.dist_el = obj.locations.getmeandist;
                obj.locations.dist_el_bi = obj.locations.getmeanbipolardist;
                %unipolar
                d_geo = obj.getdistance_geodaesic_mean();%obj.getdistance_geodaesic_mean();
                d_euclid = obj.locations.getmeandist();
                %geodisic dist should always be > than euclidean
                %but only use geodesic for CS->SPI
                mat_id_cs2spi = obj.locations.getcs2spiids_uni; 
                d_diff = d_geo-d_euclid;
                d_diff_masked = d_diff.*mat_id_cs2spi;
                [id_r,id_c]=find(d_diff_masked>0);%for these vals geo>euclid
                d_final = d_euclid;
                d_final(id_r,id_c)=d_geo(id_r,id_c);
                d_final(id_c,id_r)=d_geo(id_c,id_r);
                
                obj.locations.dist_el_geodesic = d_final;
                
                %bipolar
                d_geo_bi = obj.getdistance_geodaesic_mean_bipolar;
                %[cs_loc_bi,spi_loc_bi] = obj.locations.getbipolarlocations;
                %d_euclid_bi = squareform(pdist([squeeze(mean(cs_loc_bi,1)),squeeze(mean(spi_loc_bi,1))]'));
                d_euclid_bi = obj.locations.dist_el_bi;
                %geodisic dist should always be > than euclidean
                if size(obj.signal.bi_split.SPI,2)==9 || size(obj.signal.bi_split.SPI,2)==18 || size(obj.signal.bi_split.SPI,2)==19  
                    mat_id_cs2spi_bi = [ones(4,1); zeros(size(obj.signal.bi_split.SPI,2)+(size(d_euclid_bi,1)-4-size(obj.signal.bi_split.SPI,2)),1)] * [zeros(4,1); ones(size(obj.signal.bi_split.SPI,2)+(size(d_euclid_bi,1)-4-size(obj.signal.bi_split.SPI,2)),1)]';
                else
                    mat_id_cs2spi_bi = [ones(4,1); zeros(size(obj.signal.bi_split.SPI,2),1)] * [zeros(4,1); ones(size(obj.signal.bi_split.SPI,2),1)]';
                end
                d_diff_bi = d_geo_bi-d_euclid_bi;
                if ~isequal(size(mat_id_cs2spi_bi),size(d_diff_bi))
                    d_final_bi = d_geo_bi;
                else
                    fprintf('Not sure when this case applies. I think it has to do with CARTO, when you only have unipolar data? NEEDS A LOOKING AT!!\n')
                    waitforbuttonpress
                    d_diff_bi_masked = d_diff_bi.*mat_id_cs2spi_bi;
                    [id_r_bi,id_c_bi]=find(d_diff_bi_masked>0);
                    d_final_bi = d_euclid_bi;
                    d_final_bi(id_r_bi,id_c_bi)=d_geo_bi(id_r_bi,id_c_bi); 
                    d_final_bi(id_c_bi,id_r_bi)=d_geo_bi(id_c_bi,id_r_bi);
                end
                
                %fprintf('Check ME!!!!\n')
%                 for ii=1:size(d_geo_bi,1)
%                     for jj=1:size(d_geo_bi,1)
%                         if ii>jj
%                             d_geo_bi(ii,jj)=0; %either 0 or copy upper matrix
%                         end
%                     end
%                 end
                
                
                obj.locations.dist_el_geodesic_bi = d_final_bi;
                
            elseif isempty(obj.map) && ~isempty(obj.locations)
                obj.locations.dist_el = obj.locations.getmeandist;
                
                %add 5% of way to CS->Spi to adjust for curvature error?
                
                
                %d_euclid = obj.locations.getmeandist();
                %d_final = d_euclid;
            elseif ~isempty(obj.map) && isempty(obj.locations)
%                 inp_el_spaceing_spi = input('Input electrodespacing for spiral cath (e.g 1-4-1): ','s');
%                 el_spaceing_spi_str = regexp(inp_el_spaceing_spi,'-','split');
%                 el_spaceing_spi(1,1) = min(str2double(el_spaceing_spi_str));%short distance between elecs;
%                 el_spaceing_spi(1,2) = max(str2double(el_spaceing_spi_str));%longest
%                 
%                 inp_el_spaceing_cs = input('Input electrodespacing for cs cath (e.g 1-4-1): ','s');
%                 el_spaceing_cs_str = regexp(inp_el_spaceing_cs,'-','split');
%                 el_spaceing_cs(1,1) = min(str2double(el_spaceing_cs_str));%short distance between elecs;
%                 el_spaceing_cs(1,2) = max(str2double(el_spaceing_cs_str));
%                 
%                 obj.getdistance_spi_theo(obj,elec_spaceing);%TODO
%                 obj.getdistance_cs_theo(obj,elec_spaceing);%TODO
            elseif isempty(obj.map) && isempty(obj.locations)
                %obj.getdistance_theo(obj,elec_spaceing);
            end
        end
        
        function plotLassoCVRestitution_real(obj)
            fprintf('CARE. USES EUCLIDEAN dist only, not GEODAESIC')
            if obj.stimchan.idpostsplit-1==0
                left_of_stimchan = obj.stimchan.size_chan;
                right_of_stimchan = obj.stimchan.idpostsplit+1;
            elseif obj.stimchan.idpostsplit+1==obj.stimchan.size_chan+1
                left_of_stimchan = obj.stimchan.idpostsplit;
                right_of_stimchan = 1;
            else
                right_of_stimchan = obj.stimchan.idpostsplit + 1;
                left_of_stimchan = obj.stimchan.idpostsplit - 1;
            end
            obj.stimchan.exclude_chan = unique([obj.stimchan.exclude_chan;obj.stimchan.idpostsplit;right_of_stimchan;left_of_stimchan]);
            
            if strcmpi(obj.stimchan.cath,'CS')
                obj.stimchan.exclude_chan = obj.stimchan.idpostsplit;
            end
            %differentiate if location is given or theoretical 20mm should
            %be used
            if isempty(obj.locations) %use theoretical distance of 20mm
                fprintf('Warning: No locationdata found, using d_theo = 20mm .\n')
                d_theo=20;
                
                stim2cs_dist_all_times=1;
                stim2spi_dist_all_times = 1;
                
                stim2cs = stim2cs_dist_all_times(:,1:2:end);
                stim2spi = stim2spi_dist_all_times(:,1:end-2); %-1 because last 2 electrodes are ABL catheter & not used to detect atrial activity
                
                
                obj.atrial_segments.s2_aa_dist_fs(find(obj.atrial_segments.s2_aa_dist_fs==0)) = NaN;
                for chan=1:size(obj.atrial_segments.s2_activseg,1)
                    if chan~=obj.stimchan.idpostsplit & chan~=obj.stimchan.idpostsplit+1 & chan~=obj.stimchan.idpostsplit-1
                        figure
                        hold on
                        plot(d_theo./obj.atrial_segments.s2_aa_dist_fs(chan,:)*obj.samplerate)
                        titlestring = ['Chan:' num2str(chan)];
                        title(titlestring)
                        hold off
                    elseif chan==abs(obj.stimchan.size_chan-(floor(obj.stimchan.idpostsplit+obj.stimchan.size_chan/2)))
                        figure
                        hold on
                        plot(d_theo./obj.atrial_segments.s2_aa_dist_fs(chan,:)*obj.samplerate)
                        titlestring = ['Chan:' num2str(chan)];
                        title(titlestring)
                        hold off
                    end
                end
            else
                %Find left & right channel TODO will be taken out _> find
                %radius that is effected by stimulus and choose electrodes
                %based on that
                [stim2cs_dist_all_times, stim2spi_dist_all_times] = obj.locations.getdistances_postsplit(obj.stimchan);
                
                stim2cs = stim2cs_dist_all_times(:,1:2:end);
                if size(stim2spi_dist_all_times,2)==11 || size(stim2spi_dist_all_times,2)==21
                    stim2spi = stim2spi_dist_all_times(:,1:end-2); %-1 because last 2 electrodes are ABL catheter & not used to detect atrial activity
                else
                    stim2spi = stim2spi_dist_all_times;
                end
                
                if strcmp(obj.stimchan.cath,'CS')
                    stimchan_data = stim2cs;
                    nonstimchan_data = stim2spi;
                    stimchan_ids_of_channels = ~ismember([1:size(stim2cs,2)],obj.stimchan.idpostsplit);
                    stimchan_channels = [1:numel(stimchan_ids_of_channels)];
                    stimchan_channels = stimchan_channels(stimchan_ids_of_channels~=0);
                    %nonstimchan_channels = [1:size(stim2spi,2)];
                    all_elecs_bi = cell2mat(cellfun(@str2double,regexp(obj.elecs.bi_split.SPI,'-','split'),'UniformOutput',false));
                    all_elecs_bi(isnan(all_elecs_bi))=1; %in case a D was not removed
                    nonstimchan_channels = unique(all_elecs_bi(:,1))';
                elseif strcmp(obj.stimchan.cath,'SPI')
                    stimchan_data = stim2spi;
                    nonstimchan_data = stim2cs;
                    stimchan_ids_of_channels = ~ismember([1:size(stim2spi,2)],left_of_stimchan) .* ~ismember([1:size(stim2spi,2)],right_of_stimchan) .* ~ismember([1:size(stim2spi,2)],obj.stimchan.idpostsplit);
                    stimchan_channels = [1:numel(stimchan_ids_of_channels)];
                    stimchan_channels = stimchan_channels(stimchan_ids_of_channels~=0);
                    nonstimchan_channels = [1:size(stim2cs,2)];
                    nonstimchan_name = 'CS';
                end
                
                
                %Plotting
                obj.atrial_segments.s2_aa_dist_fs(find(obj.atrial_segments.s2_aa_dist_fs==0)) = NaN;
                obj.atrial_segments.s2_aa_dist_fs = obj.atrial_segments.s2_aa_dist_fs(:,1:obj.erp_id-1);
                
                x = mean(stimchan_data,1); %TODO for now, later every s2 time seperatly
                x_plot = obj.classcheck.s2times(1:size(obj.atrial_segments.s2_aa_dist_fs,2));
                x_plot = x_plot(1:obj.erp_id-1);
                
                stimchan_channels = [1:obj.stimchan.size_chan];
                stimchan_channels(obj.stimchan.exclude_chan) = [];
                
                clearvars name
                fig_cv_stim = figure;
                ax_cv_stim = axes;
                obj.cv_plot = nan(numel(obj.classcheck.s2times),obj.stimchan.size_chan); %init for save
                cnt=0;
                colors = jet(obj.stimchan.size_chan);
                legend_list = [];
                for chan=1:size(stimchan_channels,2)
                    
                        cnt=cnt+1;
                        %name{chan} = [obj.stimchan.cath ' ' num2str(find(x==0)) ' - ' obj.stimchan.cath ' ' num2str(stimchan_channels(chan)) ' ' 'bipolar'];
                        name_c{cnt} = [obj.stimchan.cath ' ' num2str(find(x==0)) ' - ' obj.stimchan.cath ' ' num2str(stimchan_channels(chan)) ' ' 'bipolar'];
                        hold on
                        index = stimchan_channels(chan);
                        
                        y_val = x(index)./obj.atrial_segments.s2_aa_dist_fs(index,:)*obj.samplerate;
                        y_val(obj.bad_aa_blocks) = nan;
                        
                        xs = x_plot(~isnan(y_val));
                        ys = y_val(~isnan(y_val));
                        if numel(xs)>=2
                            yi = interp1(xs, ys, x_plot, 'Linear');
                            legend_list = [legend_list 1];
                            plot(ax_cv_stim,xs,ys,'xk','LineWidth',3,'HandleVisibility','off') %to show measurementpoints
                            plot(ax_cv_stim,x_plot,yi,'color',colors(chan,:),'LineWidth',3)
                        elseif numel(xs)==0
                            yi = [];
                            legend_list = [legend_list 0];
                        else
                            yi = ys;
                            legend_list = [legend_list 1];
                            plot(ax_cv_stim,xs,ys,'xk','LineWidth',3,'HandleVisibility','off') %to show measurementpoints
                            plot(ax_cv_stim,x_plot,yi,'color',colors(chan,:),'LineWidth',3)
                        end
                        
                        %plot(ax_cv_stim,x_plot,yi,'-',xs,ys,'xk','color',colors(chan,:),'LineWidth',3)
                        %savestruct{chan} = [x_plot y_val'];
                        
                        %raw values
                        %plot(x_plot,y_val,'LineWidth',3)
                        hold off
                        %title(name)
                        %xlim([1 600]);
                        %ylim([1 1000]);
                        
                        %obj.cv_plot(1:obj.erp_id-1,stimchan_channels(chan)) = y_val(1:obj.erp_id-1)';
                   
                end
                legend(name_c(find(legend_list==1)))
                title('Stimchannel CV')
                
                
%                 clearvars name
%                 figure
%                 
%                 for chan=1:size(stimchan_channels,2)
%                     
%                     name{chan} = [obj.stimchan.cath ' ' num2str(find(x==0)) ' - ' obj.stimchan.cath ' ' num2str(stimchan_channels(chan)) ' ' 'bipolar'];
%                     hold on
%                     index = stimchan_channels(chan);
%                     
%                     y_val = squeeze(obj.atrial_segments.s2_aa_list(chan,:,2)-obj.atrial_segments.s2_aa_list(chan,:,1) );%s2_aa_dist_fs(index,:)*obj.samplerate;
%                     y_val(obj.bad_aa_blocks) = nan;
%                     
%                     xs = x_plot(~isnan(y_val));
%                     ys = y_val(~isnan(y_val));
%                     yi = interp1(xs, ys, x_plot, 'Linear');
%                     
%                     %linearly interpolated data
%                     hold on
%                     plot(xs,ys,'-xk','LineWidth',3,'HandleVisibility','off') %to show measurementpoints
%                     plot(x_plot,yi,'LineWidth',3)
%                     
%                     %savestruct{chan} = [x_plot y_val'];
%                     
%                     %raw values
%                     %plot(x_plot,y_val,'LineWidth',3)
%                     hold off
%                     %title(name)
%                     %xlim([1 600]);
%                     %ylim([1 1000]);
%                     
%                     %obj.cv_plot(1:obj.erp_id-1,stimchan_channels(chan)) = y_val(1:obj.erp_id-1)';
%                 endBES
%                 legend(name)
%                 title('Activationduration')
                
                
                
                
                
                
                
                
                %Plotting nonstimchannel
                obj.atrial_segments.nonstim_s2_aa_dist_fs(find(obj.atrial_segments.nonstim_s2_aa_dist_fs==0)) = NaN;
                obj.atrial_segments.nonstim_s2_aa_dist_fs = obj.atrial_segments.nonstim_s2_aa_dist_fs(:,1:obj.erp_id-1);
                
                x = mean(nonstimchan_data,1); %TODO for now, later every s2 time seperatly
                x_plot = obj.classcheck.s2times(1:size(obj.atrial_segments.nonstim_s2_aa_dist_fs,2));
                x_plot = x_plot(1:obj.erp_id-1);
                
                %nonstimchan_channels = [1:size(obj.signal.bi_split.(obj.stimchan.non_stimchan),2)];
                %nonstimchan_channels(obj.stimchan.exclude_chan) = []
                %exclusion of chennels TODO
                
                clearvars name
                fig_cv_nonstim = figure;
                ax_cv_nonstim = axes;
                obj.cv_plot = nan(numel(obj.classcheck.s2times),obj.stimchan.size_chan); %init for save
                cnt=0; legend_list=[];
                colors = jet(size(nonstimchan_channels,2));
                for chan=1:size(nonstimchan_channels,2)
                    cnt=cnt+1;
                    %name{chan} = [obj.stimchan.cath ' ' num2str(obj.stimchan.idpostsplit) ' - ' obj.stimchan.non_stimchan ' ' num2str(nonstimchan_channels(chan)) ' ' 'bipolar'];
                    name_c{cnt} = [obj.stimchan.cath ' ' num2str(obj.stimchan.idpostsplit) ' - ' obj.stimchan.non_stimchan ' ' num2str(nonstimchan_channels(chan)) ' ' 'bipolar']; %obj.elecs.bi_split.(obj.stimchan.cath)(obj.stimchan.idpostsplit)
                        
                    %if ~ismember(obj.stimchan.exclude_chan,chan)
                    hold on
                    index = nonstimchan_channels(chan);
                    
                    y_val = x(index)./obj.atrial_segments.nonstim_s2_aa_dist_fs(chan,:)*obj.samplerate;
                    %y_val(obj.bad_aa_blocks) = nan;
                    
                    xs = x_plot(~isnan(y_val));
                    ys = y_val(~isnan(y_val));
                    if numel(xs)>=2
                        yi = interp1(xs, ys, x_plot, 'Linear');
                        %linearly interpolated data
                        hold on
                        plot(xs,ys,'-xk','LineWidth',3,'HandleVisibility','off') %to show measurementpoints
                        plot(x_plot,yi,'color',colors(chan,:),'LineWidth',3)
                        hold off
                        %title(name)
                        %end
                        legend_list = [legend_list 1];
                    elseif numel(xs)==0
                        yi = [];
                        legend_list = [legend_list 0];
                    else
                        yi = ys;
                        %linearly interpolated data
                        hold on
                        plot(xs,ys,'-xk','LineWidth',3,'HandleVisibility','off') %to show measurementpoints
                        plot(x_plot,yi,'color',colors(chan,:),'LineWidth',3)
                        hold off
                        legend_list = [legend_list 1];
                        %title(name)
                        %end
                    end
                    
                    
                end
                legend(name_c(find(legend_list==1)))
                title('Non-Stimchannel CV')
            end
            
            
            
        end
        
        function plotLassoAMPRestitution(obj)
            
            if obj.stimchan.idpostsplit-1==0
                left_of_stimchan = obj.stimchan.size_chan;
                right_of_stimchan = obj.stimchan.idpostsplit+1;
            elseif obj.stimchan.idpostsplit+1==obj.stimchan.size_chan+1
                left_of_stimchan = obj.stimchan.idpostsplit;
                right_of_stimchan = 1;
            else
                right_of_stimchan = obj.stimchan.idpostsplit + 1;
                left_of_stimchan = obj.stimchan.idpostsplit - 1;
            end
            obj.stimchan.exclude_chan = unique([obj.stimchan.exclude_chan;obj.stimchan.idpostsplit;right_of_stimchan;left_of_stimchan]);
            
            if strcmpi(obj.stimchan.cath,'CS')
                obj.stimchan.exclude_chan = obj.stimchan.idpostsplit;
            end
            
            %Plots Amplituderestitution curve ignoring bad s2 blocks
            x_plot = obj.classcheck.s2times(1:size(obj.atrial_segments.s2_aa_dist_fs,2));
            x_plot = x_plot(1:obj.erp_id-1);
            
            stimchan_channels = [1:obj.stimchan.size_chan];
            stimchan_channels(obj.stimchan.exclude_chan) = [];
            nonstimchan_channels = cell2mat(cellfun(@str2double,regexp(obj.elecs.bi_split.(obj.stimchan.non_stimchan),'-','split'),'UniformOutput',false));
            nonstimchan_channels = nonstimchan_channels(:,1);
            %Nan vals will be interpolated linearly to make plot nicer
            fig_amp_stim = figure;
            ax_amp_stim = axes;
            colors = jet(size(stimchan_channels,2));
            cnt=0;
            legend_list = [];
            for chan=1:size(stimchan_channels,2)
                %if ~ismember(chan,exclude_chan)
                index = stimchan_channels(chan);
                    cnt=cnt+1;
                    hold on
                    vals_of_blocks = obj.signal_properties.s2_aa_p2p(index,:);
                    vals_of_blocks(obj.bad_aa_blocks) = nan;
                    vals_of_blocks = vals_of_blocks(1:obj.erp_id-1);
                    
                    xs = x_plot(~isnan(vals_of_blocks));
                    ys = vals_of_blocks(~isnan(vals_of_blocks));
                    if numel(xs)>=2
                        yi = interp1(xs, ys, x_plot, 'Linear');
                        %linearly interpolated data
                        hold on
                        plot(xs,ys,'-xk','LineWidth',3,'HandleVisibility','off') %to show measurementpoints
                        plot(x_plot,yi,'color',colors(chan,:),'LineWidth',3)
                        hold off
                        %title(name)
                        %end
                        legend_list = [legend_list 1];
                    elseif numel(xs)==0
                        yi = [];
                        legend_list = [legend_list 0];
                    else
                        yi = ys;
                        %linearly interpolated data
                        hold on
                        plot(xs,ys,'-xk','LineWidth',3,'HandleVisibility','off') %to show measurementpoints
                        plot(x_plot,yi,'color',colors(chan,:),'LineWidth',3)
                        hold off
                        legend_list = [legend_list 1];
                        %title(name)
                        %end
                    end
                    
                    %plot(x_plot,vals_of_blocks(1:obj.erp_id-1),'LineWidth',3)
                    %plot(flipud(fliplr(obj.signal_properties.s2_aa_p2p(chan,:))),'LineWidth',3)
                    %name_c{cnt} = ['Chan: ' num2str(index)];
                    name_c{cnt} = [obj.stimchan.cath ' ' obj.stimchan.idpostsplit ' - ' obj.stimchan.cath ' ' num2str(stimchan_channels(chan)) ' ' 'bipolar'];
                        
                %end
            end
            legend(name_c(find(legend_list==1)))
            hold off
            title('Stimchannel Amp')
            
            obj.amp_plot = nan(obj.erp_id-1,1);
            obj.amp_plot = vals_of_blocks(1:obj.erp_id-1)';
            
            
            
            %Plotting nonstimchannel
            %Nan vals will be interpolated linearly to make plot nicer
            fig_amp_nonstim = figure;
            ax_amp_nonstim = axes;
            legend_list = [];
            colors = jet(size(obj.signal.bi_split.(obj.stimchan.non_stimchan),2));
            for chan=1:size(obj.signal.bi_split.(obj.stimchan.non_stimchan),2)
                
                    hold on
                    vals_of_blocks = obj.signal_properties.nonstim_s2_aa_p2p(chan,:);
                    vals_of_blocks(obj.bad_aa_blocks) = nan;
                    vals_of_blocks = vals_of_blocks(1:obj.erp_id-1);
                    
                    xs = x_plot(~isnan(vals_of_blocks));
                    ys = vals_of_blocks(~isnan(vals_of_blocks));
                    if numel(xs)>=2
                        yi = interp1(xs, ys, x_plot, 'Linear');
                        %linearly interpolated data
                        hold on
                        plot(xs,ys,'-xk','LineWidth',3,'HandleVisibility','off') %to show measurementpoints
                        plot(x_plot,yi,'color',colors(chan,:),'LineWidth',3)
                        hold off
                        %title(name)
                        %end
                        legend_list = [legend_list 1];
                    elseif numel(xs)==0
                        yi = [];
                        legend_list = [legend_list 0];
                    else
                        yi = ys;
                        %linearly interpolated data
                        hold on
                        plot(xs,ys,'-xk','LineWidth',3,'HandleVisibility','off') %to show measurementpoints
                        plot(x_plot,yi,'color',colors(chan,:),'LineWidth',3)
                        hold off
                        legend_list = [legend_list 1];
                        %title(name)
                        %end
                    end
                    %plot(x_plot,vals_of_blocks(1:obj.erp_id-1),'LineWidth',3)
                    %plot(flipud(fliplr(obj.signal_properties.s2_aa_p2p(chan,:))),'LineWidth',3)
                    name_c{chan} =  [obj.stimchan.cath ' ' num2str(obj.stimchan.idpostsplit) ' - ' obj.stimchan.non_stimchan ' ' num2str(nonstimchan_channels(chan)) ' ' 'bipolar'];
                    %['Chan: ' num2str(chan)];
            end
            legend(name_c(find(legend_list==1)))
            hold off
            title('Non-Stimchannel Amp')
        
        end
        
        function dist_p2p_geodaesic = getdistance_geodaesic_mean(obj)
            %gets geodesic mean from unipolar electrodelocations if map is
            %given. Adds distance from cs to 3dmap to the pathlength
            meanloc = squeeze(mean(obj.locations.xyz,1))';
            dist_p2p_geodaesic = zeros(size(obj.locations.xyz,3));
            for i=1:size(obj.locations.xyz,3)
                for j=i+1:size(obj.locations.xyz,3)
                    [k1(i,j),d_k1(i,j)] = dsearchn(obj.map.xyz,meanloc(i,:));
                    [k2(i,j),d_k2(i,j)] = dsearchn(obj.map.xyz,meanloc(j,:));
                    [dist_p2p_geodaesic(i,j),c0{i,j},d0{i,j}] = connectPointsInMesh(obj.map.xyz,obj.map.faces,[k1(i,j); k2(i,j)]);
                    %need extrapath since cs mesh is not included
                    if strcmpi(obj.locations.cath(i),'cs') && strcmpi(obj.locations.cath(j),'cs')
                        extra_path = 0; %d_k1(i,j) + d_k2(i,j); will trvel dist along coronary sin
                    elseif strcmpi(obj.locations.cath(i),'cs') && strcmpi(obj.locations.cath(j),'spi')
                        extra_path = d_k1(i,j);
                    elseif strcmpi(obj.locations.cath(i),'spi') && strcmpi(obj.locations.cath(j),'cs')
                        extra_path = d_k2(i,j);
                    elseif strcmpi(obj.locations.cath(i),'spi') && strcmpi(obj.locations.cath(j),'spi')
                        extra_path = 0;
                    end
                    dist_p2p_geodaesic(i,j) = dist_p2p_geodaesic(i,j) + extra_path;
                end
            end
            
            %make symmetric
            for i=1:size(obj.locations.xyz,3)
                for j=i+1:size(obj.locations.xyz,3)
                    dist_p2p_geodaesic(j,i) = dist_p2p_geodaesic(i,j);
                end
            end
%             figure;
%             hold all;
%             TRtmp = triangulation(double(obj.map.faces), double(obj.map.xyz));
%             trimesh(TRtmp, 'edgecolor', [.8 .8 .8],'FaceAlpha',0.5);
%             pos = squeeze(mean(obj.locations.xyz,1));
%             plot3(pos(1,:),pos(2,:),pos(3,:),'*b','LineWidth',2.5)
%             strtid = 1;
%             endid = size(obj.locations.xyz,3)-4; %-4 because ABL
%             for i=strtid:endid 
%                 for j=i+1:endid
%                     if ~isnan(d0{i,j}(1,:))
%                         %plot3(d0{i,j}(1,:), d0{i,j}(2,:), d0{i,j}(3,:), 'ko');
%                     end
%                     if ~isnan(c0{i,j}(1,:))
%                         plot3(c0{i,j}(1,:), c0{i,j}(2,:), c0{i,j}(3,:), 'g','LineWidth',2.5);
%                     end
%                 end
%             end
           
%             legend('true control points', 'WeightedBeziere');
%             hold off;
%             axis equal

%             figure
%             labels=cellstr(num2str([1:size(meanloc,1)]'));
%             plot3(meanloc(:,1),meanloc(:,2),meanloc(:,3),'o')
%             text(meanloc(:,1),meanloc(:,2),meanloc(:,3),labels,'VerticalAlignment','bottom','HorizontalAlignment','right')

        %necessary for all times??? probably not
        end
        
        function dist_p2p_geodaesic = getdistance_geodaesic_mean_bipolar(obj)
            [cs_loc_bi,spi_loc_bi] = obj.locations.getbipolarlocations;
            meanloc = [squeeze(mean(cs_loc_bi,1)),squeeze(mean(spi_loc_bi,1))]';
            dist_p2p_geodaesic = zeros(size(meanloc,1));
            caths = [obj.leads.bi_split.CS;repmat({'SPI'},[size(meanloc,1)-size(obj.leads.bi_split.CS,1) 1])];
            for i=1:size(meanloc,1)
                for j=i+1:size(meanloc,1)
                    [k1(i,j),d_k1(i,j)] = dsearchn(obj.map.xyz,meanloc(i,:));
                    [k2(i,j),d_k2(i,j)] = dsearchn(obj.map.xyz,meanloc(j,:));
                    dist_p2p_geodaesic(i,j) = connectPointsInMesh(obj.map.xyz,obj.map.faces,[k1(i,j); k2(i,j)]);
                    if strcmpi(caths(i),'cs') && strcmpi(caths(j),'cs')
                        extra_path = d_k1(i,j) + d_k2(i,j);
                    elseif strcmpi(caths(i),'cs') && strcmpi(caths(j),'spi')
                        extra_path = d_k1(i,j);
                    elseif strcmpi(caths(i),'spi') && strcmpi(caths(j),'cs')
                        extra_path = d_k2(i,j);
                    elseif strcmpi(caths(i),'spi') && strcmpi(caths(j),'spi')
                        extra_path = 0;
                    end
                    dist_p2p_geodaesic(i,j) = dist_p2p_geodaesic(i,j) + extra_path;
                    dist_p2p_geodaesic(j,i) = dist_p2p_geodaesic(i,j); %to make the matrix symmetrical
                end
            end
            
        end
        
        
        function segment_manually(obj)
            time_list_stimcath = [];
            fprintf('Try to segment the atrial signals very cleanly with least amount of overlap possible.\n')
            
            b = input('How many Stimulationblocks? (Integer number): ');
            s1 = input('How many S1 in one block? (Integer number): ');
            s2 = input('How many S2 in one block? (Integer number): ');
            erp_id = input('First block with no atrial signal after S2? (Integer number): ');
            
            obj.signal.work = ECG_Baseline_Removal(obj.signal.in_bi,obj.samplerate,0.01,0.5);
            obj.getbipolarsignalsandleads();    %use bipolar signals for peakdetection
            %            obj.getunipolarsignalsandleads();
            obj.splitcatheters();
            obj.getstimchannel();
            obj.signal.bi_stimchan = obj.signal.bi_split.(obj.stimchan.cath);
            
            manual_segments__stimchan_list = segmentegmtool(obj.signal.bi_stimchan,obj.leads.bi_split.(obj.stimchan.cath)); %units: samples
            
            exc = input('Chose wich channels (integer id) to exclude. (Inputformat: [1 2 3]): ');
            all_chan = 1:obj.stimchan.size_chan;
            active_chan = all_chan(~ismember(all_chan,exc));
            new_num_chan = obj.stimchan.size_chan - numel(exc);
            
            fprintf('Dumping backupfile to disk.\n')
            save('dump_man_seg_raw.mat','manual_segments__stimchan_list')
            fprintf('Finished.\n')
            
            msslo = reshape(manual_segments__stimchan_list(:,1),[],erp_id-1);
            msslo_stim = repmat(reshape(msslo(1:3,:),3,1,erp_id-1),1,new_num_chan,1); %stimulus defined by 3 vals: start mid end;
            msslo = reshape(msslo(4:end,:),3,new_num_chan,erp_id-1); %each atrial signal defined by 3 vals
            %msslo(1dim:start mid end,2dim:channels;3dim s2 train stimnumber)
            
            s1s2_peak_dt = squeeze(msslo(2,:,:)-msslo_stim(2,:,:)); %rows: active_chan, columns s2 train number
            
            obj.S1_blocknum = b;
            obj.manual_segments_stimchan_list = msslo;
            obj.numof_S1_in_one_block = s1;
            obj.erp_id = erp_id;
            
        end
        
        
        function segment_manually_rough(obj)
            %segment bipolar stimulation catheter first
            %Roughly segment signal to distinguish between stimulation and
            %atrial signal. NLEO is then used to detect peaks of signals
            
            time_list_stimcath = [];
            fprintf('Try to segment the atrial signals very cleanly with least amount of overlap possible.\n')
            
            load('manual_seg_stimcan_16.mat');
            
            obj.hp_filter_signals();            %eliminates baselinewander;
            obj.getbipolarsignalsandleads();    %use bipolar signals for peakdetection
            %            obj.getunipolarsignalsandleads();
            obj.splitcatheters();
            obj.getstimchannel();
            obj.signal.bi_stimchan = obj.signal.bi_split.(obj.stimchan.cath);
            %manual_segments__stimchan_list = segmentegmtool(obj.signal.bi_stimchan,obj.leads.bi_split.(obj.stimchan.cath)); %units: samples
            
            %b = input('How many Stimulationblocks? (Integer number): ');
            %s1 = input('How many S1 in one block? (Integer number): ');
            %s2 = input('How many S2 in one block? (Integer number): ');
            %erp_id = input('First block with no atrial signal after S2? (Integer number): ');
            
            maxnum = size(manual_segments__stimchan_list,1)/2;
            msslo = zeros(maxnum,3); %manual_segments__stimchan_list ordered
            template = [ones(numel(1:s1)*2,1); ones(numel(s2)*2,1)*2]; %s1 = 1,s2 = 2,atrialsignal_s1 = 3, atrialsignal_s2 = 0
            template(2:2:end) = 3;
            template(find(template==2)+1) = 0;
            
            cnt = 0;
            cnt_t = 0;
            for i=1:2:size(manual_segments__stimchan_list,1)
                cnt = cnt+1;
                cnt_t = cnt_t + 1;
                
                msslo(cnt,1) = manual_segments__stimchan_list(i,1); %start_fs
                msslo(cnt,2) = manual_segments__stimchan_list(i+1,1); %stop_fs
                msslo(cnt,3) = template(cnt_t);
                
                if cnt_t == numel(template)
                    cnt_t = 0;
                end
            end
            msslo = round(msslo);
            
            cnt_s1 = 0;
            cnt_s2 = 0;
            cnt_aa_s1 = 0;
            cnt_aa_s2 = 0;
            for chan=1:obj.stimchan.size_chan
                for i = 1:size(msslo,1)
                    if msslo(i,3) == 1 %s1
                        cnt_s1 = cnt_s1 + 1;
                        id_s1{chan,cnt_s1} = [msslo(i,1) msslo(i,2)];
                        %nl_s1{chan,cnt_s1} = nleo(obj.signal.bi_stimchan(msslo(i,1):msslo(i,2),chan),obj.samplerate,1);
                        %[~,nl_s1_as(chan,cnt_s1)] = getActiveSegmentsFromNLEOInShortWindow(nl_s1{chan,cnt_s1},obj.samplerate,3);
                        if cnt_s1 == b*s1
                            cnt_s1 = 0;
                        end
                    elseif msslo(i,3) == 2 %s2
                        cnt_s2 = cnt_s2 + 1;
                        id_s2{chan,cnt_s2} = [msslo(i,1) msslo(i,2)];
                        nl_s2{chan,cnt_s2} = nleo(obj.signal.bi_stimchan(msslo(i,1):msslo(i,2),chan),obj.samplerate,1);
                        [~,nl_s2_as(chan,cnt_s2)] = max(nl_s2{chan,cnt_s2});
                        %[~,nl_s2_as(chan,cnt_s2)] = getActiveSegmentsFromNLEOInShortWindow(nl_s2{chan,cnt_s2},obj.samplerate,1);
                        if cnt_s2 == b*s2
                            cnt_s2 = 0;
                        end
                    elseif msslo(i,3) == 3 %atrial activity of an s1
                        cnt_aa_s1 = cnt_aa_s1 + 1;
                        id_s1_aa{chan,cnt_aa_s1} = [msslo(i,1) msslo(i,2)];
                        %nl_s1_aa{chan,cnt_aa_s1} = nleo(obj.signal.bi_stimchan(msslo(i,1):msslo(i,2),chan),obj.samplerate,1);
                        %[~,nl_s1_aa_as(chan,cnt_aa_s1)] = getActiveSegmentsFromNLEOInShortWindow(nl_s1_aa{chan,cnt_aa_s1},obj.samplerate,3);
                        if cnt_aa_s1 == b*s1
                            cnt_aa_s1 = 0;
                        end
                    elseif msslo(i,3) == 0 %atrial activity of an s2
                        cnt_aa_s2 = cnt_aa_s2 + 1;
                        id_s2_aa{chan,cnt_aa_s2} = [msslo(i,1) msslo(i,2)];
                        %nl_s2_aa{chan,cnt_aa_s2} = nleo(obj.signal.bi_stimchan(msslo(i,1):msslo(i,2),chan),obj.samplerate,1);
                        %[~,nl_s2_aa_as(chan,cnt_aa_s2)] = getActiveSegmentsFromNLEOInShortWindow(nl_s2_aa{chan,cnt_aa_s2},obj.samplerate,3);
                        if cnt_aa_s2 == b*s2
                            cnt_aa_s2 = 0;
                        end
                    end
                end
            end
            
            
            %             chan = 9;
            %             i = 2;
            %
            %             figure
            %             subplot(2,1,1)
            %             plot(obj.signal.bi_stimchan(id_s2_aa{chan,i}(1):id_s2_aa{chan,i}(2),chan))
            %             subplot(2,1,2)
            %             plot(nl_s2_aa{chan,i})
            %subplot(2,1,2)
            %plot(nl_s2_aa{chan,i})
            
            
            
            %get last s2 to reduce stimulus overlap by using matched filter
            %approach
            for chan=1:obj.stimchan.size_chan
                
                lasts2 = [min(id_s2{chan,erp_id-1}):max(id_s2{chan,erp_id-1})]; %last s2 ids
                processed_sig = obj.signal.bi_stimchan(:,chan);
                
                maxnum_of_stim = 2;
                %think about getting threshold from lonely s2 at end 10% of baseline or
                %somth like that.
                
                if chan==obj.stimchan.idpostsplit || chan==obj.stimchan.idpostsplit-1 || chan==obj.stimchan.idpostsplit+1
                    %atrialactivity{chan,i} = NaN; %s2 atrial activity filtered away stimulus
                    %nl_aa(chan,i) = {NaN};
                    %aa_loc_global(chan,i) = NaN;
                else
                    
                    for i=1:erp_id-1
                        
                        stimsig = squeeze(obj.signal.bi_stimchan(lasts2,chan)); %last s2
                        
                        idstart(chan,i) = min(id_s2{chan,i}); %start s2 stim
                        idend(chan,i) = max(id_s2_aa{chan,i}); %end of window including atrial activity
                        
                        stimAAsig{chan,i} = obj.signal.bi_stimchan(idstart(chan,i):idend(chan,i),chan);
                        [~,stimAAsig_as(chan,i)] = getActiveSegmentsFromNLEOInShortWindow(nleo(stimAAsig{chan,i},obj.samplerate,1),obj.samplerate,3);
                        
                        %TODO upsample stimulussignal before crosscorr
                        [crl{i}, crl_l{i}] = xcorr(stimAAsig{chan,i},stimsig);
                        [maxcor_val, maxcor_id]= max(crl{i});
                        x_lagg(i,1) = crl_l{i}(maxcor_id);
                        if x_lagg(i,1) < 1
                            idstart(chan,i) = idstart(chan,i)+x_lagg(i,1)-1;
                            stimAAsig{chan,i} = obj.signal.bi_stimchan(idstart(chan,i):idend(chan,i),chan);
                            [crl{i}, crl_l{i}] = xcorr(stimAAsig{chan,i},stimsig);
                            [maxcor_val, maxcor_id]= max(crl{i});
                            x_lagg(i,1) = crl_l{i}(maxcor_id);
                        end
                        
                        %calc scale to n-th maximum
                        [pks_sig, lcs_sig] = findpeaks( obj.signal.bi_stimchan( min(id_s2{chan,i}):max(id_s2{chan,i}),chan ) );
                        pks_sorted_sig = sortrows([pks_sig lcs_sig],'descend');
                        pk_choice_sig = pks_sorted_sig(maxnum_of_stim,:);
                        
                        s2sig_amp(i) = obj.signal.bi_stimchan(min(id_s2{chan,i})+pk_choice_sig(2)-1,chan);
                        
                        [pks_stim, lcs_stim] = findpeaks( stimsig );
                        pks_sorted_stim = sortrows([pks_stim lcs_stim],'descend');
                        pk_choice_stim = pks_sorted_stim(maxnum_of_stim,:);
                        
                        stimsig_amp(i) = stimsig(pk_choice_stim(2));
                        
                        scaleingfactor(i) = s2sig_amp(i)/stimsig_amp(i);
                        
                        %scale & subtract
                        zeropad_stim = zeros(size(stimAAsig{chan,i}));
                        range = ( x_lagg(i,1):(x_lagg(i,1)+size(stimsig,1)-1) ) + 1;
                        zeropad_stim( range ) = stimsig * scaleingfactor(i);
                        subtracted{chan,i} = stimAAsig{chan,i} - zeropad_stim;
                        processed_sig(idstart(chan,i):idend(chan,i)) = subtracted{chan,i};
                        
                        baseline = mean(obj.signal.bi_stimchan(min(id_s2{chan,i})-11:min(id_s2{chan,i})-1,chan)); %10 previous samples
                        ids = idstart(chan,i);
                        ide = idstart(chan,i) + stimAAsig_as{chan,i}(2)-1;
                        ide = min(id_s2_aa{chan,i});
                        processed_sig(ids:ide) = baseline;
                        
                        %new NLEO for atrial activity
                        %[~,nl_aa(chan,i)] = getActiveSegmentsFromNLEOInShortWindow(nleo(processed_sig(idstart(chan,i):idend(chan,i)),obj.samplerate,1),obj.samplerate,3);
                        
                        atrialactivity{chan,i} = processed_sig(idstart(chan,i):idend(chan,i)); %s2 atrial activity filtered away stimulus
                        [~,nl_aa(chan,i)] = getActiveSegmentsFromNLEOInShortWindow(nleo(atrialactivity{chan,i},obj.samplerate,1),obj.samplerate,1); %get activity window after filtering
                        
                        if isempty(nl_aa{chan,i}) || size(nl_aa{chan,i},1)<2
                            [~, nl_loc(chan,i)] = max(nleo(atrialactivity{chan,i},obj.samplerate,1));
                            
                        elseif size(nl_aa{chan,i},1) == 1
                            [~, nl_loc(chan,i)] = max(nleo(atrialactivity{chan,i},obj.samplerate,1));
                            %nl_loc(chan,i) = nl_loc(chan,i) -1;
                        else
                            [~, nl_loc(chan,i)] = max(nleo(atrialactivity{chan,i}(nl_aa{chan,i}(2,1):nl_aa{chan,i}(2,2)),obj.samplerate,1));
                            %nl_loc(chan,i) = nl_aa{chan,i}(2,1) + nl_loc(chan,i) -1;
                        end
                        nl_aa{chan,i} = nl_aa{chan,i} + idstart(chan,i) - 1;
                        aa_loc_global(chan,i) = nl_loc(chan,i) + idstart(chan,i);
                        nl_s2_as_global(chan,i) = nl_s2_as(chan,i)+idstart(chan,i)-1;
                        aa_s2_dist(chan,i) = abs(aa_loc_global(chan,i) - (nl_s2_as(chan,i)+idstart(chan,i)-1));
                    end
                end
                
            end
            
            
            %get unipolar max dV/dt from bipolar window using
            dvdt_uni_s2 = zeros(chan,erp_id-1);
            %filter unipolar signal
            ECG_High_Low_Filter(obj.signal.uni_split.(obj.stimchan.cath),obj.samplerate,30,270);
            obj.signal.uni_stimchan = Notch_FilterWithoutHarmonics(ECG_High_Low_Filter(obj.signal.uni_split.(obj.stimchan.cath),obj.samplerate,30,270),obj.samplerate,50,1);
            %obj.signal.uni_stimchan = obj.signal.uni_split.(obj.stimchan.cath);
            for chan=3:obj.stimchan.size_chan
                for i=1:erp_id-1
                    idstart(chan,i) = min(id_s2{chan,i});
                    idend(chan,i) = max(id_s2{chan,i});
                    sig_uni_diff_s2{chan,i} = diff(obj.signal.uni_stimchan(idstart(chan,i):idend(chan,i),chan));
                    [~,dvdt_uni_s2(chan,i)] = max(sig_uni_diff_s2{chan,i});
                    dvdt_uni_s2(chan,i) = dvdt_uni_s2(chan,i) + idstart(chan,i) - 1;
                    
                    %use maximum of NLEO as minimum for unipolar (just a hunch but this might work out best)
                    idstart(chan,i) = aa_loc_global(chan,i);
                    idend(chan,i) = max(id_s2_aa{chan,i});
                    sig_uni_diff_aa{chan,i} = diff(obj.signal.uni_stimchan(idstart(chan,i):idend(chan,i),chan));
                    [~,dvdt_uni_aa(chan,i)] = max(sig_uni_diff_aa{chan,i});
                    dvdt_uni_aa(chan,i) = dvdt_uni_aa(chan,i) + idstart(chan,i) - 1;
                    %id_s2{chan,cnt_s2};
                    
                    aa_s2_dist_uni(chan,i) = dvdt_uni_aa(chan,i) - dvdt_uni_s2(chan,i) + 1;
                end
            end
            
            %check unipolar result
            chan=7;
            figure
            hold on
            plot(obj.signal.uni_stimchan(:,chan))
            for i=1:erp_id-1
                plot(dvdt_uni_s2(chan,i),0.1,'r-*')
                plot(dvdt_uni_aa(chan,i),0.1,'g-*')
                %check how far away NLEO max is
                plot(aa_loc_global(chan,:),0.1,'b-*')
            end
            hold off
            
            %check unipolar restitution curves
            for chan=1:9
                figure
                plot(20./aa_s2_dist_uni(chan,:)*obj.samplerate)
            end
            
            
            %run NLEO on each segment in each channel to get rough equivalent to
            %automatic detection
            
            %check global vals
            chan=7;
            figure
            hold on
            plot(obj.signal.bi_stimchan(:,chan))
            plot(aa_loc_global(chan,:),0.1,'r-*')
            plot(nl_s2_as_global(chan,:),0.1,'r-*')
            
            
            chan=7;
            for i=1:erp_id-1
                figure
                subplot(3,1,1)
                plot(stimAAsig{chan,i})
                subplot(3,1,2)
                plot(atrialactivity{chan,i})
                subplot(3,1,3)
                hold on
                plot(nleo(atrialactivity{chan,i},obj.samplerate,1))
                plot(nl_loc(chan,i),max(nleo(atrialactivity{chan,i},obj.samplerate,1)),'r-*')
                hold off
            end
            
            for chan=3:9
                figure
                plot(20./aa_s2_dist(chan,:)*obj.samplerate)
            end
            
        end
        
        
        
        function plot_s2signalsegments_with_stim_single_chan(obj,chan)
            s2_list_all = obj.timeline.getS2list_chan_all;
            for i=1:size(s2_list_all,2)
                figure
                plot(obj.signal.bi_stimchan(s2_list_all(chan,i,1):s2_list_all(chan,i,2),chan))
            end
        end
        
        
    end
end

