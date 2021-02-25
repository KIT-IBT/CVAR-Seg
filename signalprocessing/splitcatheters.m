% -------------------------------------------------------
%
%    splitcatheters - based on lead names
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

function [sig_split,leads_split,elecs_split,ids] = splitcatheters(signal,leads,elecs,system)
%SPLITCATHETERS splits cathetersignals based on theire leadnames into
%different parts
%signal contains columwise channels and rowwise samples/time
%leads are ordered the same order as signalchannels
%system gives flag wich strings should be expected (for Navx e.g.
%CS,ABL,SPI...)
%
%[sig_split,leads_split] = splitcatheters(signal,leads,system)
%
%signal : columns channels, rows timesteps
%leads : rows are names
%system : 'navx' ('carto','thythmia' not implemented yet)

if strcmp(system,'navx')
    tmpnames = names();
    namelist = tmpnames.cathnames;
    %namelist = {'ABL';'CS';'SPI'}; %standard catheter names for NavX    
    for i = 1:size(namelist,1)
        ids{i} = ~cellfun('isempty',strfind(cellstr(leads),namelist{i,1}));
        sig_split.(namelist{i}) = signal(:,ids{i});
        leads_split.(namelist{i}) = leads(ids{i},:);
        elecs_split.(namelist{i}) = elecs(ids{i},:);
    end
elseif strcmp(system,'rhythmia')
   %Todo
elseif strcmp(system,'carto')
   %Todo 
end

elecs_renamed = elecs;
if sum(cellfun(@(x) ismember('-',x),elecs)) == 0 %no '-' found meaning elecs dont have form 3-4,4-5, etc..
    loc_D = strfind(elecs,'D');
    for i=1:numel(loc_D)
        if ~isempty(loc_D{i})
            elecs_renamed{i} = '1';
        end
        exist_num = str2double(elecs_renamed{i});
        calc_num = exist_num + 1;
        
        elecs_renamed{i} =  [elecs_renamed{i} '-' num2str(calc_num)];
    end
    
    
else
    %do nothing for now
end

end

