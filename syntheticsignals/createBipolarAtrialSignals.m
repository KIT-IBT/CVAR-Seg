% -------------------------------------------------------
%
%    createBipolarAtrialSignals - creates synthethic bipolar signals
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

function [pulse_final] = createBipolarAtrialSignals(sigwidth,samplerate,flag)
%CREATEBIPOLARATRIALSIGNALS gives a list of nice signals created with the
%the moving dipole function

el_dist = 3; %mm
cv = el_dist/sigwidth*1000 * 5;
%samplerate = 1000;

switch flag
    case 'updown'
        
    case 'downup'
        el_start = 500;
        el_end = el_start+el_dist;
        pulse = createBipolarAtrialSigTemplate_dipole([750;-750],[380;+110],[el_start;-el_dist/2],[el_start;el_dist/2],1,samplerate,1); % createBipolarAtrialSigTemplate_dipole([650;100],[0;-16.5],[el_start;-el_dist/2],[el_start;el_dist/2],1,samplerate,1);
%       [100;-100]*4,
        %         [stepfun,segments] = getActiveSegmentsFromNLEOInShortWindow(nleo(pulse(:,2),samplerate,1,10),samplerate,1);
%         figure
%         hold on
%         plot(pulse(:,1),pulse(:,2))
%         plot(pulse(:,1),stepfun)
        %pulse_final = [pulse(segments{1,1}(1):segments{1,1}(2),1)*1000-mean(pulse(segments{1,1}(1):segments{1,1}(2),1)*1000) pulse(segments{1,1}(1):segments{1,1}(2),2)];
        pulse = pulse - pulse(1); %set first point to 0 so no jumps when inserting to signal;
        [~,sig_cent] = max(diff(pulse(:,2)));
        sig_width_win = 2*50 /1000;
        segments{1,1} = [pulse(sig_cent,1)-sig_width_win  pulse(sig_cent,1)+sig_width_win] *samplerate;
        pulse_final = [pulse(segments{1,1}(1):segments{1,1}(2),1)*1000-mean(pulse(segments{1,1}(1):segments{1,1}(2),1)*1000) pulse(segments{1,1}(1):segments{1,1}(2),2)];
        figure
        hold on
        plot(pulse_final(:,1),pulse_final(:,2))
        %TODO center signal !!!
    case 'pass'
        
    case 'valley'
        el_start = 100;
        el_end = el_start+el_dist;
        pulse = createBipolarAtrialSigTemplate_dipole([cv;0],[0;0],[el_start;0],[el_end;0],1,samplerate,1);
        [stepfun,segments] = getActiveSegmentsFromNLEOInShortWindow(nleo(pulse(:,2),samplerate,1,0.1),samplerate,1);
        figure
        hold on
        plot(pulse(:,1),pulse(:,2))
        plot(pulse(:,1),stepfun)
        pulse_final = [pulse(segments{1,1}(1):segments{1,1}(2),1)*1000-mean(pulse(segments{1,1}(1):segments{1,1}(2),1)*1000) pulse(segments{1,1}(1):segments{1,1}(2),2)];
        figure
        hold on
        plot(pulse_final(:,1),pulse_final(:,2))
    otherwise
        
end


end

