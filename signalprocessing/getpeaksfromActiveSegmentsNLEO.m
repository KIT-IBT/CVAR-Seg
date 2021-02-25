% -------------------------------------------------------
%
%    getActiveSegmentsFromNLEOInShortWindow  - Finds active segments 
%       based on NLEO in very short windows without moving average
%
%    Ver. 0.9.0
%
%    Created:                 Tobias Oesterlein  (15.5.2013)
%    Based on work by:       Marc Aubreville    (03.7.2009)
%    Last modified:           Tobias Oesterlein  (3.9.2013)
%
%    Institute of Biomedical Engineering
%    Universitaet Karlsruhe (TH)
%
%    http://www.ibt.kit.edu
%
%    Copyright 2000-2013 - All rights reserved.
%
% ------------------------------------------------------


function [peaks_sampl,peaks_time] = getpeaksfromActiveSegmentsNLEO(signal_nleo,activeseg,samplerate)
% gets peaks from the active segment sof the NLEO.

peaks_sampl = zeros(size(activeseg,1),1);

for i=1:size(activeseg,1)
    if ~isempty(activeseg) && size(activeseg,2)~=1 
    y_max = max(signal_nleo(activeseg(i,1):activeseg(i,2)));
    x_local_seg = find(signal_nleo(activeseg(i,1):activeseg(i,2))==y_max);
    peaks_sampl(i) = x_local_seg + activeseg(i,1);
    end
end

peaks_time = peaks_sampl./samplerate;



end

