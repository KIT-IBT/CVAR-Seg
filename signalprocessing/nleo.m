% -------------------------------------------------------
%
%    nleo  - Simplified and faster version of the NLEO operator
%
%    Ver. 0.1.0 (alpha)
%
%    Created:       Marc Aubreville        (03.07.2009)
%    Last modified: Tobias Oesterlein      (21.04.2014)
%
%    Institute of Biomedical Engineering
%    Universitaet Karlsruhe (TH)
%
%    http://www.ibt.uni-karlsruhe.de
%
%    Copyright 2000-2009 - All rights reserved.
%
% ------------------------------------------------------
%
% [signalOut,stepFunc] = nleo(signal, sampleFrequency, mode, [ k ], [ lengthOfWin ])
%
% Simplified and faster version of the NLEO operator
%
%        signal:          The signal. (N samples x C channels)
%
%        sampleFrequency: Samplingrate of the signal
%
%        mode:            0 = return unfiltered value of NLEO's output
%                         1 = lowpass NLEO's output with gaussian lowpass
%                             filter
%
%        k:               k-Factor, weighting the standard deviation in 
%                         the calculation of thresholds (usually 0.1, def: 0.1)
%
%        lengthOfWin:     Length of the window used for thresholds in
%                         seconds. (AFib: 1, Aflut: 2, SR: 5, def: 1)
%                       
%
%        


% This code is based on work by Minh P. Nguyen

function [signalOut,stepFunc] = nleo(signalOrg,sampleFrequency, mode,k, lengthOfWin)

% length of padding window
paddingValue=round(2*sampleFrequency);

% retrieve numer of electrodes in signal
nSigs=size(signalOrg,2);

% padding with first and last sample
initialPadding=ones(paddingValue,nSigs).*repmat(signalOrg(1,:),paddingValue,1);
endingPadding=ones(paddingValue,nSigs).*repmat(signalOrg(end,:),paddingValue,1);
signalOrg = [initialPadding; signalOrg; endingPadding];

% prepare memory
signalOutFinal=zeros(size(signalOrg));


for nSigCounter=1:nSigs
    
    signal=signalOrg(:,nSigCounter)';

% Load FIR lowpass filter
%load Lowpass1

% create a gaussian lowpass filter
%BT = 0.12;  % 3-dB bandwidth-symbol time
%OF = 6;     % Oversampling factor (i.e., number of samples per symbol)
%NT = 4;     % 2 symbol periods to the filters peak.
%h = gaussfir(BT,NT,OF);

%TP = HdLowpass1;

%h = fdesign.lowpass('N,F3dB',2,100,sampleFrequency);
%TP = design(h,'butter');

order_prefilter=1;
order_afterfilter=2;
lowpass_frequency=sampleFrequency/8;




% WARNING: filtfilt is not causal!!!!!!!
[z,p,k]=butter(order_prefilter,2*lowpass_frequency/sampleFrequency);
sos=zp2sos(z,p,k);
Bs = sos(:,1:3); % Second Order Section numerator coeffitiens
As = sos(:,4:6); % Second Order Section denominator coeffitiens
for i=1:size(Bs,1)
	filteredsignal = filtfilt(Bs(i,:),As(i,:),signal); % recursive filtering using SOS
end

% use filter Butterworth will lead to 10 ms phase shift
% h  = fdesign.lowpass('N,F3dB',order,lowpass_frequency,sampleFrequency);
% Hd = design(h,'butter');
% filteredsignal = filter(Hd,signal);
    
    



% Filter signal with gaussian lowpass filter
%sig(1,:) = filter(TP,signal(1,:));
%sig = sig(1,ceil(length(TP.Numerator)/2):length(sig));
%sig = cat(2,sig,zeros(1,ceil(length(TP.Numerator)/2)));

%sig=filteredsignal;        % for external filtering only


% Calculate NLEO's output
nleoOut = filteredsignal(1:end-1).^2 - ([0 filteredsignal(1:end-2)] .* filteredsignal(2:end));
nleoOut(1) = nleoOut(2);
nleoOut(end+1) = nleoOut(end);     % for external filtering only

if mode == 0
    signalOutTmp = abs(nleoOut);
end
if mode == 1
    %     % lowpass NLEO's output with gaussian lowpass filter
    %     l = length(h);
    %     signalOut = conv(abs(nleoOut),h); %%% should be taken abs or not?
    %     signalOut = signalOut(floor(l/2)+1:...
    %         length(signalOut));
    %     gap = length(signalOut) - length(nleoOut);
    %     for i = 1:gap
    %         signalOut(length(nleoOut)+1) = [];
    %     end
    
    
%     % WARNING: filtfilt is not causal!!!!!!!
    lowpass_frequency=24;
    [z,p,k]=butter(order_afterfilter,2*lowpass_frequency/sampleFrequency);
    sos=zp2sos(z,p,k);
    Bs = sos(:,1:3); % Second Order Section numerator coeffitiens
    As = sos(:,4:6); % Second Order Section denominator coeffitiens
    for i=1:size(Bs,1)
        filteredsignal = filtfilt(Bs(i,:),As(i,:),nleoOut); % recursive filtering using SOS
    end

    
    % use filter Butterworth will lead to phase shift by 10 msec
%     lowpass_frequency=24;
%     h  = fdesign.lowpass('N,F3dB',order,lowpass_frequency,sampleFrequency);
%     Hd = design(h,'butter');
%     filteredsignal = filter(Hd,nleoOut);
    
    signalOutTmp = abs(filteredsignal);
    
end

% save signal result
signalOutFinal(:,nSigCounter)=signalOutTmp';

end


signalOut=signalOutFinal(paddingValue+1:end-paddingValue,:);

if (nargin < 4)
    stepFunc=0; % abort if no step desired
end
return
