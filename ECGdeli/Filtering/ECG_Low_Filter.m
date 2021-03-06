% -------------------------------------------------------
%
%    ECG_Low_Filter  - This function uses a Fourier-based algorithm to linearly lowpass filter a signal
%
%    Ver. 1.0.0
%
%    Created:           Nicolas Pilia (22.06.2020)
%    Last modified:     Claudia Nagel (29.06.2020)
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
% [filteredsignal]=ECG_Low_Filter(signal,samplerate,lowpass_frequency,varargin)
% This function uses a Fourier-based algorithm to linearly lowpass filter a signal
%
%
% The cut-off frequencies can be freely choosen for better performance.
% The user can also choose among three different filters for the low pass 
% filter. Use the third input of the function to define the type of filter. 
% Possible filters are: 'Smooth', 'Gauss', and 'Butterworth'. The default 
% type is Butterworth.
%
% Inputs:
%       signal: ecg signal to be filtered. In the multichannel case every column
%           of the signal matrix should be a channel.
%       samplerate: sample frequency in Hz of the signal 
%       lowpass_frequency: cutoff frequency of the lowpass filter 
%       varargin: type of filter 'butterworth', 'gauss', 'smooth'

%
% Outputs:
%       filtered_signal: filtered ecg signal.
%
%
% Example Usage:
%       [filteredsignal]=ECG_Low_Filter(signal,250,150)
%
% Revision history:
%
%


function [filteredsignal]=ECG_Low_Filter(signal,samplerate,lowpass_frequency,varargin)

%Check if input signal is a standing vector. In case of a lying vector the
%signal is corrected

if isvector(signal)
    if size(signal,2)>size(signal,1)
        signal=signal';
        transposeflag=true;
    else
        transposeflag=false;
    end
else
    transposeflag=false;
end
        
%In case no filter type is given. Butterworth is used.
if isempty(varargin)
    varargin='Butterworth';
end

if strcmp(varargin,'Smooth') || strcmp(varargin,'smooth') || strcmp(varargin,'s') || strcmp(varargin,'S')
    case_var=1;
elseif strcmp(varargin,'Gauss') || strcmp(varargin,'gauss') || strcmp(varargin,'G') || strcmp(varargin,'g')
    case_var=2;
elseif strcmp(varargin,'Butterworth') || strcmp(varargin,'butterworth') || strcmp(varargin,'B') || strcmp(varargin,'b')
    case_var=3;
else
    error('Filter type not recognized');
end

%class type of signal
flagsingle=0;
if ~isa(signal,'double')
    signal=double(signal);
    flagsingle=1;
end

%Number of channels
NCH=size(signal,2);

%Extend signal to avoid bordering artifacts
l=round(samplerate*10);
filteredsignal=[zeros(l,NCH);signal;zeros(l,NCH)];
for i=1:NCH
    extension=wextend(1,'sp0',signal(:,i),l,'r');
    filteredsignal(:,i)=wextend(1,'sp0',extension,l,'l');
end

%Cut-Off frecuencies for the filters.
%lowpass_frequency=40;
if lowpass_frequency>samplerate/2
    disp('Warning: Lowpass frequency above Nyquist frequency. Nyquist frequency is chosen instead.')
    lowpass_frequency=floor(samplerate/2-1);
end

% lowpass filter
switch case_var
    case 1 %Smoothing filter
        nw=round(samplerate/lowpass_frequency);
        for i=1:NCH
            filteredsignal(:,i)=smooth(filteredsignal(:,i),nw);
        end
    case 2 %Gaussian filter
        sigmaf=lowpass_frequency;
        sigma=samplerate/(2*pi*sigmaf);
        length_gauss=2*round(4*sigma)+1;
        h=fspecial('gaussian',[length_gauss,1],sigma);
        for i=1:NCH
            filteredsignal(:,i)=conv(filteredsignal(:,i),h,'same');
        end
    case 3 %Butterworth filter
        order=3;
        [z,p,k]=butter(order,2*lowpass_frequency/samplerate);
        sos=zp2sos(z,p,k);
        Bs = sos(:,1:3); % Second Order Section numerator coeffitiens
        As = sos(:,4:6); % Second Order Section denominator coeffitiens
        for j=1:NCH
            for i=1:size(Bs,1)
                filteredsignal(:,j) = filtfilt(Bs(i,:),As(i,:),filteredsignal(:,j)); % recursive filtering using SOS
            end
        end
end

%Remove extention of signal
filteredsignal=filteredsignal(l+1:end-l,:);

%Constant offset removal
filteredsignal=Isoline_Correction(filteredsignal);
            
if flagsingle
    filteredsignal=single(filteredsignal);
end

%Transpose signal if necesary
if transposeflag
    filteredsignal=filteredsignal';
end

