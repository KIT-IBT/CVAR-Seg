% -------------------------------------------------------
%
%    location  - class for storing geometric data
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
% name          -name of the class
% samplerate    -4th dim, in case of catheters it's time
% xyz           -xyz in columns for all catheters (datpoints x 3 x cath)
% el            -electrodes {'D' '1' '2' ...}
% cath          -cathetertype {'CS', 'SPI', 'ABL'}
% time          -explicit time. If not given, it is calculated from
%                samplerate. If given it overrides samplerate
%
% Outputs:
%           location class
%
%
% Example Usage:
%    leads = {'CS';'CS';'CS';'CS';'SPI';'SPI';'SPI';'SPI';'SPI';'SPI';'SPI';'SPI';'SPI'};
%    elecs = {'1-2';'3-4';'5-6';'7-8';'1-2';'2-3';'3-4';'4-5';'5-6';'6-7';'7-8';'8-9';'9-10'};
%    el_pos = rand(18,3); %18 electrodes
%    samplerate = 500;
%
%    xyzin = permute(repmat(el_pos, [1 1 100]),[3 2 1]);
%       loc_in = location('benchmark',samplerate,xyzin,elecs,leads);
%
% Revision history:

classdef location < handle
    
    properties
        name
        samplerate
        xyz
        el
        cath
        time
        dist_el                 %Matrix of euclidean distances between all electrodes, orderded aftert 'cath'
        dist_el_bi              %row & colum showing how the electrodes where ordered when creating dist_el
        dist_el_geodesic
        dist_el_geodesic_bi
    end
    
    methods
        function obj = location(namein,sampleratein,xyzin,elin,cathin,timein)
            switch nargin
                case 0
                    
                case 1
                    if ischar(namein)
                        obj.name = namein;
                    else
                        fprintf('location class: Name should be string.\n');
                    end
                case 2
                    if ischar(namein)
                        obj.name = namein;
                    else
                        fprintf('location class: Name should be string.\n');
                    end
                    if isnumeric(sampleratein) && size(sampleratein,1) == 1 && size(sampleratein,2) == 1
                        obj.samplerate = sampleratein;
                    end
                    obj.time = (0:size(xyzin,1)-1)/sampleratein;
                case 3
                    if ischar(namein)
                        obj.name = namein;
                    else
                        fprintf('location class: Name should be string.\n');
                    end
                    if isnumeric(sampleratein) && size(sampleratein,1) == 1 && size(sampleratein,2) == 1
                        obj.samplerate = sampleratein;
                    end
                    if size(xyzin,1) > 0 && size(xyzin,2) == 3 && size(xyzin,3) > 0 && size(xyzin,3) > size(xyzin,1) %Last condition could cause problems with small traces
                        obj.xyz = xyzin;
                    else 
                        maxind = find(max(size(xyzin))); %Finds largest dimension, wich is probably datapoints
                        dimind = find(size(xyzin)==3);   %Finds the xyz points dimension
                        cathind = 6-maxind-dimind; %6 because all id's together 1+2+3=6 so if we subract the two other ones we get the third
                        
                        obj.xyz = permute(xyzin,[maxind dimind cathind]);
                    end
                    obj.time = (0:size(xyzin,1)-1)/sampleratein;
                case 4
                    if ischar(namein)
                        obj.name = namein;
                    else
                        fprintf('location class: Name should be string.\n');
                    end
                    if isnumeric(sampleratein) && size(sampleratein,1) == 1 && size(sampleratein,2) == 1
                        obj.samplerate = sampleratein;
                    end
                    if size(xyzin,1) > 0 && size(xyzin,2) == 3 && size(xyzin,3) > 0 && size(xyzin,3) > size(xyzin,1) %Last condition could cause problems with small traces
                        obj.xyz = xyzin;
                    else 
                        maxind = find(max(size(xyzin))); %Finds largest dimension, wich is probably datapoints
                        dimind = find(size(xyzin)==3);   %Finds the xyz points dimension
                        cathind = 6-maxind-dimind; %6 because all id's together 1+2+3=6 so if we subract the two other ones we get the third
                        
                        obj.xyz = permute(xyzin,[maxind dimind cathind]);
                    end
                    if size(elin,1) > size(elin,2)
                        obj.el = elin;
                    else
                        obj.el = elin';
                    end
                    obj.time = (0:size(xyzin,1)-1)/sampleratein;
                case 5
                    if ischar(namein)
                        obj.name = namein;
                    else
                        fprintf('location class: Name should be string.\n');
                    end
                    if isnumeric(sampleratein) && size(sampleratein,1) == 1 && size(sampleratein,2) == 1
                        obj.samplerate = sampleratein;
                    end
                    if size(xyzin,1) > 0 && size(xyzin,2) == 3 && size(xyzin,3) > 0 && size(xyzin,3) > size(xyzin,1) %Last condition could cause problems with small traces
                        obj.xyz = xyzin;
                    else 
                        maxind = find(max(size(xyzin))); %Finds largest dimension, wich is probably datapoints
                        dimind = find(size(xyzin)==3);   %Finds the xyz points dimension
                        cathind = 6-maxind-dimind; %6 because all id's together 1+2+3=6 so if we subract the two other ones we get the third
                        
                        obj.xyz = permute(xyzin,[maxind dimind cathind]);
                    end
                    if size(elin,1) > size(elin,2)
                        obj.el = elin;
                    else
                        obj.el = elin';
                    end
                    if size(cathin,1) > size(cathin,2)
                        obj.cath = cathin;
                    else
                        obj.cath = cathin';
                    end
                    obj.time = (0:size(xyzin,1)-1)/sampleratein;
                case 6
                    if ischar(namein)
                        obj.name = namein;
                    else
                        fprintf('location class: Name should be string.\n');
                    end
                    if isnumeric(sampleratein) && size(sampleratein,1) == 1 && size(sampleratein,2) == 1
                        obj.samplerate = sampleratein;
                    end
                    if size(xyzin,1) > 0 && size(xyzin,2) == 3 && size(xyzin,3) > 0 && size(xyzin,3) > size(xyzin,1) %Last condition could cause problems with small traces
                        obj.xyz = xyzin;
                    else 
                        maxind = find(max(size(xyzin))); %Finds largest dimension, wich is probably datapoints
                        dimind = find(size(xyzin)==3);   %Finds the xyz points dimension
                        cathind = 6-maxind-dimind; %6 because all id's together 1+2+3=6 so if we subract the two other ones we get the third
                        
                        obj.xyz = permute(xyzin,[maxind dimind cathind]);
                    end
                    if size(elin,1) > size(elin,2) && ndims(elin)==2
                        obj.el = elin;
                    else
                        obj.el = elin';
                    end
                    if size(cathin,1) > size(cathin,2) && ndims(cathin)==2
                        obj.cath = cathin;
                    else
                        obj.cath = cathin';
                    end
                    if size(timein,1) > size(timein,2) && ndims(timein)==2
                        obj.time = timein;
                    else
                        obj.time = timein';
                    end
            end
        end
        
        function distmatrix = getdistances(obj)
            %gets strait distances between all elecetrodes for each timestep
            %outputformat is: cath x cath x time
            distmatrix = nan(size(obj.xyz,1),size(obj.xyz,3),size(obj.xyz,3));
            
            %namematrix = nan(size(obj.xyz,1),size(obj.xyz,3),size(obj.xyz,3));
            %cathel = cellfun(@(x,y) [x y],obj.cath,obj.el,'UniformOutput',false);
            
            %Matrix aus namen erzeugen irgendwie
            
            for i = 1:size(obj.xyz,3) %runs through catheters
                distmatrix(:,i,:) = squeeze(sqrt((obj.xyz(:,1,i)-obj.xyz(:,1,:)).^2+(obj.xyz(:,2,i)-obj.xyz(:,2,:)).^2+(obj.xyz(:,3,i)-obj.xyz(:,3,:)).^2));
                %namematrix(:,i,:) = [obj.cath];
            end
            distmatrix = permute(distmatrix,[2 3 1]);
        end
        
        function calc_dist_el(obj)
            obj.dist_el = obj.getmeandist;
        end
        
        function calc_dist_bi(obj)
           obj.dist_el_bi = obj.getmeanbipolardist;
        end
        
        function [cs, spi] = getdistances_postsplit(obj,stimchan)
            idlist_cs = find(cellfun(@(x) strcmp(x,'CS'),obj.cath)==1);
            size_cs = numel(idlist_cs);
            idlist_spi = find(cellfun(@(x) strcmp(x,'SPI'),obj.cath)==1);
            size_spi = numel(idlist_spi); 
            
            cnt=1;
            for i=1:size(obj.cath,1)
                if strcmp(obj.cath{i},stimchan.cath)
                    
                    if cnt==stimchan.el1
                        stimchan_id_1 = i;
                        stimchan_id_2 = i+1;
                    end
                    cnt = cnt+1;
                end
            end
            
            %get only cs % only SPI
            %cs
            mean_pos_stim(:,1) = mean([obj.xyz(:,1,stimchan_id_1),obj.xyz(:,1,stimchan_id_2)],2);
            mean_pos_stim(:,2) = mean([obj.xyz(:,2,stimchan_id_1),obj.xyz(:,2,stimchan_id_2)],2);
            mean_pos_stim(:,3) = mean([obj.xyz(:,3,stimchan_id_1),obj.xyz(:,3,stimchan_id_2)],2);
            
            for i = 1:numel(idlist_cs)-1
                mean_pos_chan(:,1) = mean([obj.xyz(:,1,idlist_cs(i)),obj.xyz(:,1,idlist_cs(i+1))],2);
                mean_pos_chan(:,2) = mean([obj.xyz(:,2,idlist_cs(i)),obj.xyz(:,2,idlist_cs(i+1))],2);
                mean_pos_chan(:,3) = mean([obj.xyz(:,3,idlist_cs(i)),obj.xyz(:,3,idlist_cs(i+1))],2);
                cs(:,i) = squeeze(sqrt((mean_pos_chan(:,1)-mean_pos_stim(:,1)).^2+(mean_pos_chan(:,2)-mean_pos_stim(:,2)).^2+(mean_pos_chan(:,3)-mean_pos_stim(:,3)).^2));
            end
            for i = 1:numel(idlist_spi)-1
                mean_pos_chan(:,1) = mean([obj.xyz(:,1,idlist_spi(i)),obj.xyz(:,1,idlist_spi(i+1))],2);
                mean_pos_chan(:,2) = mean([obj.xyz(:,2,idlist_spi(i)),obj.xyz(:,2,idlist_spi(i+1))],2);
                mean_pos_chan(:,3) = mean([obj.xyz(:,3,idlist_spi(i)),obj.xyz(:,3,idlist_spi(i+1))],2);
                spi(:,i) = squeeze(sqrt((mean_pos_chan(:,1)-mean_pos_stim(:,1)).^2+(mean_pos_chan(:,2)-mean_pos_stim(:,2)).^2+(mean_pos_chan(:,3)-mean_pos_stim(:,3)).^2));
            end
          
        end
        
        function reduced_location_matrix = get_filtered_location_of_cath(obj,cathname)
            %extraxts locationmatrix of only the catheterelectrodes
            %specified by cathname -> e.g 'cs' gets times x cs_elecs x xyz matrix
            %cathname should be 'cs' or 'spi'
            idlist = find(cellfun(@(x) strcmpi(x,cathname),obj.cath)==1);
            reduced_location_matrix = obj.xyz(:,:,idlist);
        end
        
        function meandistmatrix = getmeandist(obj)
            meandistmatrix = mean(obj.getdistances,3); %mean over times
        end
        
        function [cs_loc_bip,spi_loc_bip] = getbipolarlocations(obj)
            
            cs_id = strcmpi(obj.cath,'cs');
            spi_id = strcmpi(obj.cath,'spi');
            
            cs_xyz = obj.xyz(:,:,cs_id);
            spi_xyz = obj.xyz(:,:,spi_id);
            
            cnt=0;
            for chan=1:2:size(cs_xyz,3)
                cnt=cnt+1;
                for i=1:size(cs_xyz,1)
                    cs_loc_bip(i,:,cnt) = mean([cs_xyz(i,:,chan);cs_xyz(i,:,chan+1)],1);
                end
            end
            cnt=0;
            for chan=1:size(spi_xyz,3)-1
                cnt=cnt+1;
                for i=1:size(spi_xyz,1)
                    spi_loc_bip(i,:,cnt) = mean([spi_xyz(i,:,chan);spi_xyz(i,:,chan+1)]);
                end
            end
            
        end
        
        function [idmatrix,idmatrix_reduced] = getcs2spiids_uni(obj)
            cs_id = strcmpi(obj.cath,'cs');
            spi_id = strcmpi(obj.cath,'spi');
            idmatrix = logical(cs_id * spi_id');
            
            if sum(cs_id)==8 %Matrix should be ordered CS then Spi
                if sum(spi_id) == 12 || sum(spi_id) == 22
                    spi_id(find(spi_id==1,2,'last'))=0;%delete last 2 elecs
                    idmatrix_reduced = logical(cs_id * spi_id');
                else
                    idmatrix_reduced = idmatrix;
                end
            else
                fprintf('TODO: Implement locationclass fuction\n')
            end
        end
        
        function [idmatrix,idmatrix_reduced] = getcs2spiids_bi(obj)
            fprintf('Not sure about this function\n')
            cs_id = [1;1;1;1];
            spi_id = strcmpi(obj.cath,'spi');
            if sum(spi_id) == 12 || sum(spi_id) == 22
                spi_id = [zeros(sum(cs_id),1); ones(sum(spi_id),1)];
                idmatrix = logical(cs_id * spi_id');
                spi_id(find(spi_id==1,3,'last'))=0;%for unipolar 12 elecs: -2 since last 2 elecs useless and -1 to convert to having bipolar leads
                idmatrix_reduced = logical(cs_id * spi_id');
            else
                idmatrix_reduced = idmatrix;
            end
        end
        
        function distmat_bi = getbipolardist(obj)
            [cs_loc_bip,spi_loc_bip] = obj.getbipolarlocations;
            comb_loc = cat(3,cs_loc_bip,spi_loc_bip);
            for i = 1:size(cs_loc_bip,3)+size(spi_loc_bip,3) %runs through catheters
                distmatrix(:,i,:) = squeeze(sqrt((comb_loc(:,1,i)-comb_loc(:,1,:)).^2+(comb_loc(:,2,i)-comb_loc(:,2,:)).^2+(comb_loc(:,3,i)-comb_loc(:,3,:)).^2));
                %namematrix(:,i,:) = [obj.cath];
            end
            distmat_bi = permute(distmatrix,[2 3 1]);
        end
        
        function distmat_bi = getmeanbipolardist(obj)
            distmat_bi = mean(obj.getbipolardist,3);
        end
        
        
    end
end

