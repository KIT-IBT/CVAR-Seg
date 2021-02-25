% -------------------------------------------------------
%
%    timeline  - class for saving stimulation segments of multichannel
%    measurements
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
% timeline(segmentarray,weights)
% class storing timesegments and corrsponding id data. 
%
% 
% Inputs:
% timeline(segmentarray,weights)
%
% Outputs:
%           timeline class
%
%
% Example Usage:
%   timeline_all = timeline();
%   weightvector = rand(1,10);
%   thresholdvector = [ones(1,5)];
%   timeline_all.initialize_weights(weightvector);      %as of now weights are given as ids referring to this array
%   timeline_all.initialize_thresholds(thresholdvector);
%   timeline_all.initialize_stimchanid(1);
%   timeline_all.addsegments(1,10,5,1,2,'Methodname')
%   timeline_all.createsegmentids();
%   timeline_all.createsupersegments(); % 2 is weightid for detection in another channel, second is threshold id
%   timeline_all.createsupersegmentids();
%   timeline_all.updateweight();
%
% Revision history:


classdef timeline < handle
    
    
    properties
        segments            
        supersegments
        numofsegments
        numofsupersegments
        
        threshold
        signallength_fs
        stimchanid
        
        distance
        std
        
        %weights
        weights
        weightstring        %comment what weight is for e.g NLEO is first weight,
        weightthreshold
        
    end
    
    properties (Access=private)
        
    end
    
    methods
        function obj = timeline(segmentarray,weights)
            %Constructs
            if nargin~=0
                p = inputParser;
                addRequired(p,'segmentarray',@(x) (isvector(x) && isa(x,'timesegment')) || isvector(x));
                addRequired(p,'weights',@(x) isvector(x) && isnumeric(x) && max(size(weights)) == max(size(segmentarray)));
                addRequired(p,'stimchanid',@(x) (isnumeric(x) && size(x,1) == 1 && size(x,2) == 1));
                try
                    p.parse( segmentarray,weights,stimchanid );
                catch MyError
                    rethrow(MyError);
                end
                
                if size(p.Results.segmentarray,1) > size(p.Results.segmentarray,2)
                    obj.segments = p.Results.segmentarray;
                else
                    obj.segments = p.Results.segmentarray';
                end
                obj.numofsegments = max(size(segments));
                obj.weights = p.Results.weights;
            else
                obj.numofsegments = 0;
                obj.weights = 0;
            end
            obj.numofsupersegments = 0;
            
        end
        
        function initialize_weights(obj,weightarray)
            if ~isempty(weightarray) && isnumeric(weightarray) && isvector(weightarray)
                if size(weightarray,2)>size(weightarray,1)
                    obj.weights = weightarray;
                else
                    obj.weights = weightarray';
                end
            else
                fprintf('wrong inputformat of weightarray\n')
            end
        end
        
        function initialize_thresholds(obj,thresholdarray)
            if ~isempty(thresholdarray) && isnumeric(thresholdarray) && isvector(thresholdarray)
                if size(thresholdarray,2)>size(thresholdarray,1)
                    obj.threshold = thresholdarray;
                else
                    obj.threshold = thresholdarray';
                end
            else
                fprintf('wrong inputformat of weightarray\n')
            end
        end
        
        function initialize_stimchanid(obj,stimchanid)
            if ~isempty(stimchanid) && isnumeric(stimchanid) && size(stimchanid,1)==1 && size(stimchanid,2)==1
                obj.stimchanid = stimchanid;
            else
                fprintf('Wrong stimchanid input\n')
            end
        end
        
        
        %add segments
        function addsinglesegmentbypeakandsiglength(obj,peak_fs,siglength,weightidofsegment,weightidifinotherchan)
            if size(peak_fs,1) == 1 && size(peak_fs,2) == 1 && size(siglength,1) == 1 && size(siglength,2) == 1
                obj.update_numofsegments;
                %check if existing segment is hit
                new_seg = timesegment(peak_fs-siglength/2,peak_fs+siglength/2,peak_fs,'weight',obj.weights(weightidofsegment));
                segment_hit_id = obj.checkifsegmentfitsexistingsegments(new_seg);
                
                if isempty(segment_hit_id)
                    %adds segment to timeline
                    obj.segments{obj.numofsegments+1} = new_seg;
                else
                    obj.segments{segment_hit_id}.weight = obj.segments{segment_hit_id}.weight + obj.weights(weightidifinotherchan);
                end
                
                %check
                
            else
                fprintf('Is not a single peak but array. doing nothing.\n')
                %call arraysegmentreadin function later on
            end
        end
        
        function addsupersegementcenter(obj,id)
            %check if overlap with existing segments
            obj.update_numofsupersegments;
            
            overlapid = obj.check_overlapwithallsupersegments(id);
            
            if isnumeric(id) && sum(size(id))==2 && isempty(overlapid)
                %obj.check_superseg_add(id)
                obj.supersegments{obj.numofsupersegments+1,1} = timesegment(obj.segments{id}.min_fs,obj.segments{id}.max_fs,obj.segments{id}.peaks);
            else
                fprintf('Supersegments overlap, ignoring supersegment.\n');
                %add supersegment to existing supersegment
            end
            if ~isempty(obj.supersegments{obj.numofsupersegments+1,1}.weight)
                obj.supersegments{obj.numofsupersegments+1,1}.weight = [];
            end
            %obj.supersegments{obj.numofsupersegments+1,1}.addcontent(timesegment(obj.segments{id}.min_fs,obj.segments{id}.max_fs,obj.segments{id}.peaks)); %adding channel itself as content as well
            obj.update_numofsupersegments;
        end
        
        function createsegmentids(obj)
            obj.update_numofsegments
            for i=1:obj.numofsegments
                obj.segments{i}.segment_id = i;
            end
        end
        
        function createsupersegmentids(obj)
            obj.update_numofsupersegments
            %obj.sortsupersegments
            for i=1:obj.numofsupersegments
                obj.supersegments{i}.segment_id = i;
            end
        end
        
        function sortsupersegments(obj)
            obj.update_numofsupersegments
            for i=1:obj.numofsupersegments
                peaklist(i) = obj.supersegments{i}.peaks(1);
            end
            [~,ind] = sort(peaklist);
            for i=1:obj.numofsupersegments
                sorted_list{i,1} = obj.supersegments{ind(i)};
            end
            obj.supersegments = sorted_list;
        end
        
        function addsegmentweightbyid(obj,id,weight)
            if ~isnumeric(id) || size(id,1) ~= 1 || size(id,2) ~= 1 || id<=0
                fprintf('id has wrong inputformat.\n')
                return
            end
            if ~isnumeric(weight) || size(weight,1) ~= 1 || size(weight,2) ~= 1
                fprintf('id has wrong inputformat.\n')
                return
            end
            
            
            if id<=obj.numofsegments
                obj.segments{id}.weight = weight;
            else
                fprintf('id is out of bounds.\n')
            end
            
        end
        
        function segment_hit_id = checkifsegmentfitsexistingsegments(obj,segment)
            %return id id segment exists and nothing if not
            %edges are part of segment
            if isa(segment,'timesegment')
                obj.update_numofsegments;
                for i = 1:obj.numofsegments
                    isinsegment(i) = obj.segments{i}.min_fs <= segment.peaks && obj.segments{i}.max_fs >= segment.peaks;
                end
                segments_hit_id = find(isinsegment==1);
                if max(size(isinsegment(isinsegment==1))) > 1
                    fprintf('More than one segment was hit. Not possible.\n')
                    return
                else
                    segment_hit_id = segments_hit_id;
                end
            else
                fprintf('Is not a segment. doing nothing.\n')
                %call apropriate function later on
            end
            
        end
        
        function checkifpeakfitsexistingsegment(obj,peak)
            %return id id segment exists and nothing if not
            %edges are part of segment
            obj.update_numofsegments;
            for i = 1:obj.numofsegments
                isinsegment(i) = obj.segments{i}.min_fs <= peak && obj.segments{i}.max_fs >= peak;
            end
            segments_hit_id = find(isinsegment==1);
            if max(size(isinsegment(isinsegment==1))) > 1
                fprintf('Multiple segments share same peaks...strange.\n')
                return
            end
        end
        
        function seg_ids = checkifpeakfitsexistingsupersegment(obj,peak)
            %return id id segment exists and nothing if not
            %edges are part of segment
            obj.update_numofsupersegments;
            for i = 1:obj.numofsupersegments
                isinsegment(i) = obj.supersegments{i}.min_fs <= peak && obj.supersegments{i}.max_fs >= peak;
            end
            segments_hit_id = find(isinsegment==1);
            if max(size(isinsegment(isinsegment==1))) > 1
                fprintf('Not possible.\n')
                seg_ids = [];
            elseif isempty(segments_hit_id)
                seg_ids = [];
            else
                seg_ids = segments_hit_id;
            end
        end
        
        function addsegments(obj,min_fs,max_fs,peaks,weightid,chanid_postsplit,description)
            obj.update_numofsegments;
            if isempty(peaks)
                numofpeaks=0;
            else
                numofpeaks = max(size(peaks));
            end
            for i=1:numofpeaks
                obj.segments{obj.numofsegments+i,1} = timesegment(min_fs(i),max_fs(i),peaks(i),'weight',obj.weights(weightid),'chanid_postsplit',chanid_postsplit);
                if ~isempty(description) && ischar(description)
                    obj.segments{obj.numofsegments+i,1}.history{1} = description;
                end
            end
            obj.check_if_min_negative;
        end
        
        function addsegmentsbypeakandsiglength(obj,peaks,guess_siglength_stim,weightid,chanid_postsplit)
            obj.update_numofsegments;
            %adds segment arrays to timeline
            numofpeaks = max(size(peaks));
            
            for i=1:numofpeaks
                obj.segments{obj.numofsegments+i,1} = timesegment(peaks(i)-guess_siglength_stim/2,peaks(i)+guess_siglength_stim/2,peaks(i),'weight',obj.weights(weightid),'chanid_postsplit',chanid_postsplit);
            end
            %check if input overlaps with existing segments
            %obj.check_for_overlapping_segments
            obj.check_if_min_negative;
        end
        
                
        %updatefunctions
        function update_numofsegments(obj)
            obj.numofsegments = size(obj.segments,1);
        end
        
        function update_numofsupersegments(obj)
            obj.numofsupersegments = size(obj.supersegments,1);
        end
        
        function updateweight(obj)
            obj.update_numofsegments
            obj.update_numofsupersegments
            currentsegweight=[];
            for i=1:obj.numofsupersegments
                for j=1:size(obj.supersegments{i}.content,1)
                    %fprintf('%d,%d\n',i,j)
                    currentsegweight(j) = obj.supersegments{i}.content{j}.weight;
                end
                if isempty(currentsegweight)
                    obj.supersegments{i}.weight = 0;
                else
                    obj.supersegments{i}.weight = sum(currentsegweight);
                end
                currentsegweight = 0;
            end
        end
        
        
        %checkfunctions
        function check_if_min_negative(obj)
            %there should be no negative values in timeline
            for i=1:obj.numofsegments
                if obj.segments{i}.min_fs < 0
                    fprintf('Negative segment found.\n')
                    obj.segments{i}.min_fs = 0;
                end
            end
        end
        
        function check_for_overlapping_segments(obj)
            obj.update_numofsegments;
            
            cnt = 0;
            for i = 1:obj.numofsegments             %i = small number -> max
                for j = i+1:obj.numofsegments       %j = larger -> min
                    cnt = cnt + 1;
                    overlap_matrix(cnt,1) = i;
                    overlap_matrix(cnt,2) = j;
                    overlap_matrix(cnt,3) = obj.segments{i}.max_fs < obj.segments{j}.min_fs;
                end
            end
            
            overlap_ids = find(overlap_matrix(:,3)~=1,1);
            
            if ~isempty(overlap_ids)
                fprintf('Overlap of segments found.\n')
                %Decide what to do with these segments once this case
                %actually appears
                return
            end
            
        end
        
        function stimchanids = getstimchansegmentids(obj,stimchanid_postsplit)
            obj.update_numofsegments
            cnt = 0;
            for i=1:obj.numofsegments
                if obj.segments{i}.channelid_postsplit == stimchanid_postsplit
                    cnt = cnt+1;
                    stimchanids(cnt) = i; 
                end
            end
        end
        
        function overlap = checkifsegmentsoverlap(obj,seg1,seg2,threshold)
            if abs(mean(seg1.peaks) - mean(seg2.peaks)) <= threshold
                overlap = true;
            else
                overlap = false;
            end
        end
        
        function [segment_overlap,segment_overlap_id_peak,segment_overlap_id_min,segment_overlap_id_max] = check_overlapwithallsupersegments(obj,segmentid)
            obj.update_numofsupersegments
            cnt_ignored = 0;
            segment_overlap_id_max = [];
            segment_overlap_id_min = [];
            segment_overlap_id_peak = [];
            segment_overlap_peak_id_array = [];
            segment_overlap_min_id_array = [];
            segment_overlap_max_id_array = [];
            segment_overlap = [];
            
            if obj.numofsupersegments > 0
                cnt=0;
                cnt_max=0;
                cnt_min=0;
                for i=1:obj.numofsupersegments
                    
                    if obj.segments{segmentid}.peaks <= obj.supersegments{i}.max_fs && obj.segments{segmentid}.peaks >= obj.supersegments{i}.min_fs
                        cnt=cnt+1;
                        segment_overlap_peak_id_array(cnt) = i;
                    elseif obj.segments{segmentid}.min_fs >= obj.supersegments{i}.min_fs &&  obj.segments{segmentid}.min_fs <= obj.supersegments{i}.max_fs %mininsegment
                        cnt_min=cnt_min+1;
                        segment_overlap_min_id_array(cnt_min) = i;
                    elseif obj.segments{segmentid}.max_fs >= obj.supersegments{i}.min_fs &&  obj.segments{segmentid}.min_fs <= obj.supersegments{i}.max_fs%maxinsegment
                         cnt_max=cnt_max+1;
                        segment_overlap_max_id_array(cnt_max) = i;
                    end
                    
                end
            end
            
            seg_ovlap_peak = segment_overlap_peak_id_array(segment_overlap_peak_id_array~=0);
            seg_ovlap_min = segment_overlap_min_id_array(segment_overlap_min_id_array~=0);
            seg_ovlap_max = segment_overlap_max_id_array(segment_overlap_max_id_array~=0);
            peak_seg_size = numel(seg_ovlap_peak);
            min_seg_size = numel(seg_ovlap_min);
            max_seg_size = numel(seg_ovlap_max);
            
            
            if peak_seg_size==1 || min_seg_size==1 || max_seg_size==1 %if one peak was found -> return val
                segment_overlap_id_peak = seg_ovlap_peak;
                segment_overlap_id_min = seg_ovlap_min;
                segment_overlap_id_max = seg_ovlap_max;
                
                if ~isempty(segment_overlap_id_peak) || ~isempty(segment_overlap_id_min) || ~isempty(segment_overlap_id_max)
                    segment_overlap = 1;
                else
                    segment_overlap = [];
                end
                
            elseif peak_seg_size>1 || min_seg_size>1 || max_seg_size>1
                fprintf('More than one segment was hit. Ignoring for now.\n')
                cnt_ignored = cnt_ignored + 1;
            else
                segment_overlap_id_peak = [];
                segment_overlap_id_min = [];
                segment_overlap_id_max = [];
                segment_overlap = [];
            end
            
        end
        
        function segment_overlap_id = check_overlapwithallsupersegments_old(obj,segmentid)
            obj.update_numofsupersegments
            
            segment_overlap_id_array = [];
            if obj.numofsupersegments > 0
                cnt=0;
                for i=1:obj.numofsupersegments
                    
                    if obj.segments{segmentid}.peaks <= obj.supersegments{i}.max_fs && obj.segments{segmentid}.peaks >= obj.supersegments{i}.min_fs
                        cnt=cnt+1;
                        segment_overlap_id_array(cnt) = i;
%                     else
%                         cnt=cnt+1;
%                         segment_overlap_id_array(cnt) = [];
                    end
                end
            end
            
            seg_ovlap = segment_overlap_id_array(segment_overlap_id_array~=0);
            if (size(seg_ovlap,1)==1 && size(seg_ovlap,2)==1) || isempty(seg_ovlap)
                segment_overlap_id = seg_ovlap;
            else
                fprintf('More than one segment was hit.\n')
               
             end
            
        end
        
%         %wrapperfunctions
        function createsupersegments(obj)
            obj.update_numofsegments;
            
            for i=1:obj.numofsegments
                
                [overlap,overlap_peak,overlap_min,overlap_max] = obj.check_overlapwithallsupersegments(obj.segments{i}.segment_id);
                %i
                if overlap
                    overlap_logic = zeros(1,3);
                    if ~isempty(overlap_peak)
                        overlap_logic(1) = 1;
                    end
                    if ~isempty(overlap_min)
                        overlap_logic(2) = 1;
                    end
                    if ~isempty(overlap_max)
                        overlap_logic(3) = 1;
                    end
                    
                    overlap_array = [overlap_peak,overlap_min,overlap_max];
                    switch numel(overlap_array)
                        case 1
                            overlapid = overlap_array;
                        case 2
                            if isequal(overlap_array(1),overlap_array(2))
                                overlapid = overlap_array(1);
                            else
                                fprintf('Peakposition is highest priority. Adding to segment.\n')
                                overlapid = overlap_array(1);
                            end
                        case 3
                            if isequal(overlap_array(1),overlap_array(2),overlap_array(3))
                                overlapid = overlap_array(1);
                            else
                                fprintf('Case not implemented either.\n')
                                %Todo for now add peak as most prominent
                                overlapid = overlap_array(1);
                            end
                        otherwise
                            fprintf('Code should not reach this point. Probably empty overlapid\n')
                            return
                    end
                    
                    obj.add2supersegment(overlapid,obj.segments{i});
                else
                    obj.addsupersegment(obj.segments{i});
                end

            end
            
        end
        
        function add2supersegment(obj,overlapid,segment)
            obj.update_numofsupersegments
            if isempty(obj.supersegments{overlapid}.content)
                obj.supersegments{overlapid}.content{1,1} = segment;
            else
                obj.supersegments{overlapid}.content{end+1,1} = segment;
            end

        end
        
        function addsupersegment(obj,segment)
            obj.update_numofsupersegments
            if isempty(segment.content)
                obj.supersegments{obj.numofsupersegments+1,1} = timesegment(segment.min_fs,segment.max_fs,segment.peaks,'weight',segment.weight);
                obj.supersegments{obj.numofsupersegments+1,1}.addcontent(segment);
            else %Should never enter this case
                obj.supersegments{obj.numofsupersegments+1,1} = timesegment(segment.min_fs,segment.max_fs,segment.peaks,'weight',segment.weight);
                for i=1:max(size(segment.content))
                    obj.supersegments{obj.numofsupersegments+1,1}.addcontent(segment.content{i});
                end
            end
            %delete weighting since this will be done in an extrafunction
            %after completion
            if ~isempty(obj.supersegments{obj.numofsupersegments+1,1}.weight)
                obj.supersegments{obj.numofsupersegments+1,1}.weight=[];
            end
            
        end
        
        function segment = getsegbyid(obj,id)
            obj.update_numofsegments
            for i=1:obj.numofsegments
                if obj.segments{i}.segment_id==id
                    segment = obj.segments{i};
                end
            end
        end
        
        function addweight2supersegmentbyid(obj,id,weight)
            if ~isnumeric(id) || size(id,1) ~= 1 || size(id,2) ~= 1 || id<=0
                fprintf('id has wrong inputformat.\n')
                return
            end
            if ~isnumeric(weight) || size(weight,1) ~= 1 || size(weight,2) ~= 1
                fprintf('id has wrong inputformat.\n')
                return
            end
            
            
            if id<=obj.numofsupersegments
                obj.supersegments{id}.weight = obj.supersegments{id}.weight + weight;
            else
                fprintf('id is out of bounds.\n')
            end
        end
        
        function addweight2supersegbypeak(obj,weight,peak)
            for i = 1:size(peak,1)
                id = obj.checkifpeakfitsexistingsupersegment;
                fprintf('NEINEINEIN')
            end
        end
        
        
        
        
        %plotting functions
        function plotsupersegments(obj,signal,weightthreshold)
            %plots segments against whatever signal is inputtet
            %signalinput should be rows - time and columns - channels
            figure;
            hold on
            obj.update_numofsupersegments
            y_max = 0;
            %plot signal
            for i=1:size(signal,2)
                plot(signal(:,i),'b')
                if max(abs(signal(:,i)))>y_max
                    y_max=max(abs(signal(:,i)));
                end
            end
            %plot segments
            segmentcount=0;
            for i=1:obj.numofsupersegments
                if obj.supersegments{i}.weight >= weightthreshold
                    patch([obj.supersegments{i}.min_fs,obj.supersegments{i}.min_fs,obj.supersegments{i}.max_fs,obj.supersegments{i}.max_fs],[-ceil(y_max),ceil(y_max),ceil(y_max),-ceil(y_max)],'red','FaceColor', [255 150 0]/255, 'FaceAlpha', .5, 'EdgeColor', 'none')
                    segmentcount=segmentcount+1;
                end
            end
            box on
            fprintf('%.0f segments found with %.0f weighting.\n',segmentcount,weightthreshold)
        end
        
        function plotsegments(obj,signal)
            %plots segments against whatever signal is inputtet
            %signalinput should be rows - time and columns - channels
            figure;
            hold on
            obj.update_numofsegments
            y_max = 0;
            %plot signal
            for i=1:size(signal,2)
                plot(signal(:,i),'b')
                if max(abs(signal(:,i)))>y_max
                    y_max=max(abs(signal(:,i)));
                end
            end
            %plot segments
            segmentcount=0;
            for i=1:obj.numofsegments
                
                    patch([obj.segments{i}.min_fs,obj.segments{i}.min_fs,obj.segments{i}.max_fs,obj.segments{i}.max_fs],[-ceil(y_max),ceil(y_max),ceil(y_max),-ceil(y_max)],'red','FaceColor', [255 150 0]/255, 'FaceAlpha', .5, 'EdgeColor', 'none')
                    segmentcount=segmentcount+1;
                
            end
            %fprintf('%.0f segments found with %.0f weighting.\n',segmentcount,weightthreshold)
        end
        
        function f1 = plotsupersegmentweights(obj,signal)
            f1 = figure;
            hold on
            yyaxis left
            ylim([0 Inf])
            obj.update_numofsupersegments
            y_max = 0;
            %plot signal
            for i=1:size(signal,2)
                plot(signal(:,i),'-b');
                if max(abs(signal(:,i)))>y_max
                    y_max=max(abs(signal(:,i)));
                end
            end
            yyaxis right
            %ylim([-max(obj.getuniqueweights) Inf])
            %plot weights
            segmentcount=0;
            for i=1:obj.numofsupersegments
                plot(obj.supersegments{i}.peaks(1,1),obj.supersegments{i}.weight,'-*','LineWidth',3)
                segmentcount=segmentcount+1;
            end
            fprintf('%.0f segments found.\n',segmentcount)
        end
        
        function plots1s2(obj,signal)
            figure;
            hold on
            obj.update_numofsupersegments
            y_max = 0;
            %plot signal
            for i=1:size(signal,2)
                plot(signal(:,i),'b')
                if max(abs(signal(:,i)))>y_max
                    y_max=max(abs(signal(:,i)));
                end
            end
            %plot s1s2
            segmentcounts1=0;
            segmentcounts2=0;
            for i=1:obj.numofsupersegments
                if strcmp(obj.supersegments{i}.segtype,'S1')
                    patch([obj.supersegments{i}.min_fs,obj.supersegments{i}.min_fs,obj.supersegments{i}.max_fs,obj.supersegments{i}.max_fs],[-ceil(y_max),ceil(y_max),ceil(y_max),-ceil(y_max)],'blue','FaceColor', [255 150 0]/255, 'FaceAlpha', .5, 'EdgeColor', 'none')
                    segmentcounts1=segmentcounts1+1;
                elseif strcmp(obj.supersegments{i}.segtype,'S2')
                    patch([obj.supersegments{i}.min_fs,obj.supersegments{i}.min_fs,obj.supersegments{i}.max_fs,obj.supersegments{i}.max_fs],[-ceil(y_max),ceil(y_max),ceil(y_max),-ceil(y_max)],'red','FaceColor', [220 20 60]/255, 'FaceAlpha', .5, 'EdgeColor', 'none')
                    segmentcounts2=segmentcounts2+1;
                end
            end
%             for i=1:obj.numofsupersegments
%                 if strcmp(obj.supersegments{i}.segtype,'S1')
%                     patch([obj.supersegments{i}.min_fs,obj.supersegments{i}.min_fs,obj.supersegments{i}.max_fs,obj.supersegments{i}.max_fs],[-ceil(y_max),ceil(y_max),ceil(y_max),-ceil(y_max)],'blue','FaceColor', [255 150 0]/255, 'FaceAlpha', .5, 'EdgeColor', 'none')
%                     segmentcounts1=segmentcounts1+1;
%                 elseif strcmp(obj.supersegments{i}.segtype,'S2')
%                     for j=1:size(obj.supersegments{i}.content,1)
%                         allmin(j) = obj.supersegments{i}.content{j}.min_fs;
%                         allmax(j) = obj.supersegments{i}.content{j}.max_fs;
%                     end
%                     patch([max(allmin),max(allmin),min(allmax),min(allmax)],[-ceil(y_max),ceil(y_max),ceil(y_max),-ceil(y_max)],'red','FaceColor', [220 20 60]/255, 'FaceAlpha', .5, 'EdgeColor', 'none')
%                     segmentcounts1=segmentcounts2+1;
%                 end
%             end
            fprintf('%.0i s1 found.\n',segmentcounts1)
            fprintf('%.0i s2 found.\n',segmentcounts2)
        end
        
        function plotwithpeaks(obj,signal,peaks)
            figure;
            hold on
            %plot signal
            y_max = 0;
            for i=1:size(signal,2)
                plot(signal(:,i),'b')
                if max(abs(signal(:,i)))>y_max
                    y_max=max(abs(signal(:,i)));
                end
            end
            %plot peaks
            for i=1:max(size(peaks))
                plot(peaks(i),1,'x','LineWidth',2)
            end
        end
        
        function plottimeline(obj,signal,timelsegments)
            
            hold on
            y_max = 0;
            %plot signal
            for i=1:size(signal,2)
                figure
                plot(signal(:,i),'b')
                if max(abs(signal(:,i)))>y_max
                    y_max=max(abs(signal(:,i)));
                end
            end
            %plot segments
            segmentcount=0;
            for i=1:max(size(timelsegments))
                if ~isempty(timelsegments{i})
                    patch([timelsegments{i}.min_fs,timelsegments{i}.min_fs,timelsegments{i}.max_fs,timelsegments{i}.max_fs],[-ceil(y_max),ceil(y_max),ceil(y_max),-ceil(y_max)],'red','FaceColor', [255 150 0]/255, 'FaceAlpha', .5, 'EdgeColor', 'none')
                    segmentcount=segmentcount+1;
                end
            end
            
        end
        
        
        
        %get functions
        function timeline_filt = gettimesegs_filtered_by_weight(obj,weighthreshold)
            obj.update_numofsupersegments
            for i=1:obj.numofsupersegments
                if obj.supersegments{i}.weight <= weighthreshold
                    timeline_filt{i} = {};
                else
                    timeline_filt{i} = obj.supersegments{i};
                end
            end
        end
        
        function idlist = get_segids_by_weight(obj,weight)
            obj.update_numofsupersegments
            cnt=0;
            for i=1:obj.numofsupersegments
                if obj.supersegments{i}.weight == weight
                    cnt = cnt+1;
                    idlist(cnt) = obj.supersegments{i}.segment_id;
                end
            end
        end
        
        function idlist = get_segids_in_weight_window(obj,min_weight,max_weight)
            obj.update_numofsupersegments
            cnt=0;
            idlist = [];
            for i=1:obj.numofsupersegments
                if obj.supersegments{i}.weight > min_weight && obj.supersegments{i}.weight < max_weight 
                    cnt = cnt+1;
                    idlist(cnt) = obj.supersegments{i}.segment_id;
                end
            end
        end
        
        function supersegids = getsupersegids_smaller_than_weight(obj,weighthreshold)
            %gets ids smaller than weight
            obj.update_numofsupersegments
            supersegids = [];
            cnt=0;
            for i=1:obj.numofsupersegments
                if obj.supersegments{i}.weight < weighthreshold
                    cnt=cnt+1;
                    supersegids(cnt) = i;
                end
            end
        end
        
        function supersegids = getsupersegids_larger_than_weight(obj,weighthreshold)
            %gets ids larger or equal than weight
            obj.update_numofsupersegments
            supersegids = [];
            cnt=0;
            for i=1:obj.numofsupersegments
                if obj.supersegments{i}.weight > weighthreshold
                    cnt=cnt+1;
                    supersegids(cnt) = i;
                end
            end
        end
        
        function uniqueweights = getuniqueweights(obj)
            obj.update_numofsupersegments;
            for i = 1:obj.numofsupersegments
                weightlist(i) = obj.supersegments{i}.weight;
            end
            uniqueweights = unique(weightlist);
        end
        
        function redmat = createreduceddistmat(obj,supersegids)
            %creates reduced matrix keeping ids of original timeline
            endsize = max(size(supersegids));
            cnt=0;
            for i=1:endsize
                for j=i+1:endsize
                    cnt=cnt+1;
                    redmat(cnt,1) = obj.supersegments{supersegids(i)}.segment_id;
                    redmat(cnt,2) = obj.supersegments{supersegids(j)}.segment_id;
                    redmat(cnt,3) = abs(obj.supersegments{supersegids(i)}.peaks(1)-obj.supersegments{supersegids(j)}.peaks(1));
                end
            end
            
        end
        
        function mean_siglength = getmeansiglength(obj)
            obj.update_numofsupersegments;
            for i=1:obj.numofsupersegments
                m(i) = abs(obj.supersegments{i}.max_fs - obj.supersegments{i}.min_fs);
            end
            mean_siglength = mean(m);
        end
        
        function median_siglength = getmediansiglength(obj)
            obj.update_numofsupersegments;
            for i=1:obj.numofsupersegments
                m(i) = abs(obj.supersegments{i}.max_fs - obj.supersegments{i}.min_fs);
            end
            median_siglength = median(m);
        end
        
        function dist = getdistbetweensupersegids(id1,id2)
            dist = abs(mean(obj.timeline.supersegments{id1}.peaks)-mean(obj.timeline.supersegments{id2}.peaks));
        end
        
        function [dist,stddiv,bad_dist_id] = getS2S1dist(obj)
            obj.update_numofsupersegments;
            
            cnt=0;
            for i=1:obj.numofsupersegments
                if strcmp(obj.supersegments{i}.segtype,'S1')
                    cnt=cnt+1;
                    list(cnt,1) = 1;
                    list(cnt,2) = i;
                elseif strcmp(obj.supersegments{i}.segtype,'S2')
                    cnt=cnt+1;
                    list(cnt,1) = 2;
                    list(cnt,2) = i;
                %else
                %    cnt=cnt+1;
                %    list(cnt) = 0;
                end
            end
            
            cnt_dist = 0;
            for i=1:cnt
                if i~=cnt %exclude if last val is 2
                    if list(i,1)==2 && list(i+1,1)~=2
                        cnt_dist = cnt_dist + 1;
                        distances(cnt_dist) = abs(obj.supersegments{i+1}.peaks(1)-obj.supersegments{i}.peaks(1));
                        idlist_dist(cnt_dist) = i;
                    end
                end
            end
            med_dist = median(distances);
            stand_div = 0.1 * med_dist; %10percent deviation
            
            dist_reduced = distances(distances<med_dist+stand_div & distances>med_dist-stand_div);
            
            missing_dist_id = find(~ismember(distances,dist_reduced));
            bad_dist_id = idlist_dist(missing_dist_id);
            
            if abs(max(size(distances))-max(size(dist_reduced))) > ceil(0.5*max(size(distances))) %if more than 50% of values were discarded
                fprintf('S2S1 distance cant be found. Distinguishment between S1S2 is not good. Bad signal?\n')
                dist=NaN;
            else
                dist = mean(dist_reduced);
                stddiv = std(dist_reduced); 
            end
            
            
        end
        
        function list = getS1list_all(obj)
            obj.update_numofsupersegments
                cnt=0;
                for i=1:obj.numofsupersegments
                    if strcmp(obj.supersegments{i}.segtype,'S1')
                        cnt=cnt+1;
                        list(cnt,1) = obj.supersegments{i}.min_fs;
                        list(cnt,2) = obj.supersegments{i}.max_fs;
                        list(cnt,3) = obj.supersegments{i}.peaks;
                    end
                end
        end
        
        function list = getS2list_all(obj)
            obj.update_numofsupersegments
                cnt=0;
                for i=1:obj.numofsupersegments
                    if strcmp(obj.supersegments{i}.segtype,'S2')
                        cnt=cnt+1;
                        list(cnt,1) = obj.supersegments{i}.min_fs;
                        list(cnt,2) = obj.supersegments{i}.max_fs;
                        list(cnt,3) = obj.supersegments{i}.peaks;
                    end
                end
        end
        
        function list = getS2list_chan_all(obj)
            obj.update_numofsupersegments
            cnt=0;
            for i=1:obj.numofsupersegments
                if strcmp(obj.supersegments{i}.segtype,'S2')
                    cnt=cnt+1;
                    for j=1:numel(obj.supersegments{i}.content)
                        id = obj.supersegments{i}.content{j}.channelid_postsplit;
                        list(id,cnt,1) = obj.supersegments{i}.content{j}.min_fs;
                        list(id,cnt,2) = obj.supersegments{i}.content{j}.max_fs;
                        list(id,cnt,3) = obj.supersegments{i}.content{j}.peaks;
                        %list(id,cnt,4) = i;
                    end
                end
            end
        end
        
        function list = getS1list_chan_all(obj)
            obj.update_numofsupersegments
            cnt=0;
            for i=1:obj.numofsupersegments
                if strcmp(obj.supersegments{i}.segtype,'S1')
                    cnt=cnt+1;
                    for j=1:numel(obj.supersegments{i}.content)
                        id = obj.supersegments{i}.content{j}.channelid_postsplit;
                        list(id,cnt,1) = obj.supersegments{i}.content{j}.min_fs;
                        list(id,cnt,2) = obj.supersegments{i}.content{j}.max_fs;
                        list(id,cnt,3) = obj.supersegments{i}.content{j}.peaks;
                    end
                end
            end
            fprintf('This function gives matrixentries with zeros and has to be checked!. AVOID IT.\n')
            %[c1,c2,c3]=find(list==0);
            %list(c1,c2,c3)=[];
        end
        
        function segments = getsupersegbystring(obj,identifier)
            if ~ischar(identifier)
                fprintf('Identifier must be string or char.\n')
                segments = NaN;
                return
            else
                obj.update_numofsupersegments
                cnt=0;
                for i=1:obj.numofsupersegments
                    if strcmp(obj.supersegments{i}.segtype,identifier)
                        cnt=cnt+1;
                        segments(cnt,1) = obj.supersegments{i};
                    end
                end
            end
        end
        
        function setsupersegbystring(obj,identifier,timesegments)
            if ~ischar(identifier)
                fprintf('Identifier must be string or char.\n')
                return
            else
                obj.update_numofsupersegments
                cnt=0;
                for i=1:obj.numofsupersegments
                    if strcmp(obj.supersegments{i}.segtype,identifier)
                        cnt=cnt+1;
                        for j=1:numel(timesegments)
                            if timesegments(j).segment_id == obj.supersegments{i}.segment_id
                                obj.supersegments{i} = timesegments(j);
                            end
                        end
                    end
                end
            end
        end
        
        function locationtable = getstimuluslocationtable(obj,identifier)
            %identifiers is array of strings to identify different peaks
            %depending on namingscheme: here: 'S1','S2'
            if ~ischar(identifier)
                fprintf('Identifier must be string or char.\n')
                locationtable = NaN;
                return
            else
                obj.update_numofsupersegments
                cnt=0;
                for i=1:obj.numofsupersegments
                    if strcmp(obj.supersegments{i}.segtype,identifier)
                        cnt=cnt+1;
                        locationtable(cnt,1) = obj.supersegments{i}.min_fs;
                        locationtable(cnt,2) = obj.supersegments{i}.max_fs;
                        locationtable(cnt,3) = mean(obj.supersegments{i}.peaks);
                    end
                end
            end
        end
        
        function locationtable = getlasts1locationtable(obj,numof_stim_in_block,identifier)
            locationtable = 0;
            s1counter = 0;
            for i=1:obj.numofsupersegments
                if strcmp(obj.supersegments{i}.segtype,identifier)
                    s1counter = s1counter + 1;
                    if s1counter==numof_stim_in_block
                        s1counter = 0;  %reset
                        row_idx = size(locationtable,1)+1;
                        locationtable(row_idx,1) = obj.supersegments{i}.min_fs;
                        locationtable(row_idx,2) = obj.supersegments{i}.max_fs;
                        locationtable(row_idx,3) = mean(obj.supersegments{i}.peaks);
                    end
                end
            end
            locationtable = locationtable(2:end,:);
        end
        
        function clearweights(obj)
            obj.update_numofsupersegments
            for i=1:obj.numofsupersegments
                obj.supersegments{i}.weight = [];
            end
        end
        
        
    end
end

