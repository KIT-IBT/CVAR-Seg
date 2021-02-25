% -------------------------------------------------------
%
%    getstimchanfromsig
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

function [stimchan] = getstimchanfromsig(signal_in,leads,elecs,type)
%GETSTIMCHANFROMSIG detects stimchannel from signal by getting channel with
%max mean amplitude and second highest meanamplitude and checks if
%neighbours
%
%signal_in should have form: rows - samples, columns - channels,
%leads is a row array of cells with one entry per cell, eg. [{'SPI'};{'SPI'};{'SPI'};{'CS'};{'CS'}]
%el is a row array of cells with one entry per cell, eg. [{'1-2'}{'2-3'}{'3-4'}] for bipolar or [{'1'};{'2'}] for unipolar
%type is a string: 'bipolar or 'unipolar'
%important: signals & leads & electrodes should be sorted the same


%get first maximum

%remove Baselinewander so maximum can be deduced)
chan_mean = mean(abs(signal_in),1); %using mean so that single noisy peaks don't destroy search
[max_val,maxid]= max(chan_mean);
stimchannels(1,1) = maxid;




if strcmp(type,'unipolar')
    %get second maximum in neighbourhood
    chan_mean(maxid) = NaN;
    chan_mean_surrounding = chan_mean(maxid-1:maxid+1);
    [~,max_channel_id2] = max(chan_mean_surrounding);
    stimchannels(2,1) = maxid-2 + max_channel_id2;
    
    %get channel from leads
    tmpnames = names();
    namelist = tmpnames.cathnames;
    stimchannel = leads{maxid};
    if ~isempty(find(stimchannel=='D'))         %Change Distal 'D' to 1
        stimchannel(stimchannel=='D')='1';
    end
    
    for i = 1:size(namelist,1)
        ids{i} = strfind(lower(stimchannel),lower(namelist{i,1})); %lower() makes it case insensitive
        if ~isempty(ids{i})
            id = i;
        end
    end
    if sum(~cellfun('isempty',ids)) ~= 1
        fprintf('getstimchanfromsig: your namingscheme contains multiple channelnames, please change it or the script.\n')
        return
    end
    stimchan.cath = namelist{id};
    
    %get postsplitid from know stimchannel
    for i = 1:size(leads,1)
        id_postid_cell{i} = logical(strfind(lower(leads{i,1}),lower(stimchan.cath)));
        id_postid(i) = ~isempty(id_postid_cell{i});
    end
    id_postid1 = find(id_postid~=0,1,'first');
    
    %unipolar result
    stimchan.el1 = min(stimchannels)-id_postid1+1;
    stimchan.el2 = max(stimchannels)-id_postid1+1;
    stimchan.idpresplit = min(stimchannels(2,1));                   %presplit is always id of first electrode
    stimchan.idpostsplit = stimchan.el2/2;
    stimchan.size_chan = size(id_postid1);
    stimchan.size_chan = abs(id_postid2-id_postid1)+1;
    

 elseif strcmp(type,'bipolar')

    %get channel from leads
    tmpnames = names();
    namelist = tmpnames.cathnames;
    stimchannel = leads{maxid};
    for i = 1:size(namelist,1)
        ids{i} = strfind(lower(stimchannel),lower(namelist{i,1})); %lower() makes it case insensitive
        if ~isempty(ids{i})
            id = i;
        end
    end
    if sum(~cellfun('isempty',ids)) ~= 1
        fprintf('getstimchanfromsig: your namingscheme contains multiple channelnames, please change it or the script.\n')
        return
    end
    stimchan.cath = namelist{id};

    %find out at what position the stimchan is to calc new ids
    for i = 1:size(leads,1)
        id_postid_cell{i} = logical(strfind(lower(leads{i,1}),lower(stimchan.cath)));
        id_postid(i) = ~isempty(id_postid_cell{i});
    end
    id_postid1 = find(id_postid~=0,1,'first');                      %first id of stimulation catheter found
    id_postid2 = find(id_postid~=0,1,'last');
    stimchansize = size(id_postid(id_postid==true),2);
    %extract second electrode from leadname
    
    replace = cellfun(@(x) strrep(x,'D','1') ,elecs,'UniformOutput',false); %Change Distal 'D' to 1
    %Is there a 'P' electrode for proxymal???
    
    
    %split electrode numbers
    split_elecs = regexp(replace,'\-','split');
    split_elecs = [cellfun(@(x) str2double(x{1,1}) ,split_elecs,'UniformOutput',true), cellfun(@(x) str2double(x{1,2}) ,split_elecs,'UniformOutput',true)];
    
    stimchan.el1 = split_elecs(maxid,1);
    stimchan.el2 = split_elecs(maxid,2);
    stimchan.idpostsplit = maxid-id_postid1+1;
    stimchan.idpresplit = maxid;
    stimchan.size_chan = abs(id_postid2-id_postid1)+1;%sum(cellfun(@(x) strcmpi(x,'CS'),leads));

end

