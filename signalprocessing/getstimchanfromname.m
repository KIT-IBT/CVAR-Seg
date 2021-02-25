% -------------------------------------------------------
%
%    getstimchanfromname 
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


function [stim_chan] = getstimchanfromname(name)
%GETSTIMCHAN gets Stimulationchannel from filename
%needs specifif naming scheme: TODO
%type: 'unipolar','bipolar'
stim_chan.cath = [];
stim_chan.el1 = [];
stim_chan.el2 = [];
stim_chan.idpostsplit = [];

tmpnames = names();
namelist = tmpnames.cathnames;

for i = 1:size(namelist,1)
    ids{i} = strfind(lower(name),lower(namelist{i,1})); %lower() makes it case insensitive
    if ~isempty(ids{i})
        id = i;
    end
end
if sum(~cellfun('isempty',ids)) ~= 1
    fprintf('getstimchanfromname: your namingscheme contains multiple channelnames, please change it or the script.\n')
    stim_chan.cath='ERROR';
    return
end

stim_chan.cath = namelist{id};


all_delimiters = strfind(lower(name),lower('_')); 
next_delim_after_stim_chan = all_delimiters(find(all_delimiters>ids{id},1,'first'))-1; %-1 since '_' should not be included

tmp = regexp(name(ids{id}:next_delim_after_stim_chan),'[A-Z]','split');
%delete everything before and after first nonempty find
tmp_idx = find(~cellfun(@isempty,tmp)~=0,1,'first');%find(~isempty(tmp),1,'first')
comment = tmp{1,tmp_idx};
numeric_ids = regexp(comment,'\d');
%numerics = tmp(1,tmp_idx);

switch size(numeric_ids,2)
    case 1
        fprintf('getstimchanfromname: Something went wrong with the naming scheme.\n')
    case 2
        stim_chan.el1 = str2double(comment(min([numeric_ids(1) numeric_ids(2)] )));
        stim_chan.el2 = str2double(comment(max([numeric_ids(1) numeric_ids(2)] )));
    case 3
        stim_chan.el1 = str2double(comment(numeric_ids(1)));
        stim_chan.el2 = str2double(comment(numeric_ids(2:3)));
    case 4
        stim_chan.el1 = str2double(comment(numeric_ids(1:2)));
        stim_chan.el2 = str2double(comment(numeric_ids(3:4)));
end

%different spaceings for CS (1-2,3-4,5-6) & Spi (1-2,2-3,3-4,4-5,...)
if strcmp(stim_chan.cath,'CS')
    stim_chan.idpostsplit = round(stim_chan.el1/2);
elseif strcmp(stim_chan.cath,'SPI')
    stim_chan.idpostsplit = stim_chan.el1;
else strcmp(stim_chan.cath,'ABL')
    fprintf('Not implemented yet.\n')
end
 
% if strcmp(type,'unipolar')
%     stim_chan.idpost = [stim_chan.el1 stim_chan.el2];
% elseif strcmp(type,'bipolar')
%     stim_chan.idpost = size(leads,1)-(size(leads,1)-stim_chan.el2+1);
% end

%get stimulus id for the block that its in instead of global id
% if strcmp(stim_chan.cath,'CS')                                   %depending where stimulus is for bipolar signal difference if stim or atrial sig
%     stim_chan.idpostsplit = ceil(stim_chan.el2/2);
% elseif strcmp(stim_chan.cath,'SPI')
%     stim_chan.idpostsplit = size(leads,1)-(size(leads,1)-stim_chan.el2+1);
% elseif strcmp(stim_chan.cath,'ABL')
%     %not implemented yet
% end

end

