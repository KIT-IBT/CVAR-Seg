% -------------------------------------------------------
%
%    addnoise_noisetemp - adds noisetemplate to the signal 
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

function [Y_time_new] = addnoise_noisetemp(signal,samplerate,noisesample,samplerate_noise,noisepwr,method,p)

% signal given in row : y x 1 
% samplerate is integer
% noisesample is either single noise sequence in timedomain or array; rows
% are function values, columns number of noisesegments
% noisepwr to be added in frequency domain;
% p is for plot
% method: check code for now; explanation following
%
% Usage: addnoise_noisetemp(zeros(2000,1),1000,rand(800,1),5,'rep',1);

%getting both signals to the same samplerate
samplerate_difference_factor = round(samplerate/samplerate_noise);
noisesample = resample(noisesample,samplerate_difference_factor,1);

Lns = length(noisesample);
Ls = length(signal);
t = (0:1/samplerate:(length(signal)-1)/samplerate);
timewin_perc_of_noisetemplate = 0.90;
switch method
    case 'rep'
        noisesample_rep = repmat(noisesample(:,1),ceil(Ls/Lns),1);
        noisesample_new = noisesample_rep;
    case 'rep_multi'
        noisesample_rep_multi = noisesample(:);
        if length(noisesample_rep_multi)<length(signal)
            noisesample_rep_multi = repmat(noisesample_rep_multi,ceil(length(signal)/length(noisesample_rep_multi)),1); %not tested
        end
        noisesample_new = noisesample_rep_multi;
    case 'randsample'
        noisesample_randsample = repmat(noisesample,ceil(Ls/Lns),1);
        timewin = ceil(timewin_perc_of_noisetemplate*Lns); %90 percent of noisewin
        samps_val = ceil(Ls/timewin);
        noisesample_new = nan(samps_val*timewin,1);
        for samps = 1:samps_val
            if samps==1
                noisesample_new(1:samps*timewin,1) = datasample(noisesample_randsample,timewin);
            else
                noisesample_new((samps-1)*timewin+1:samps*timewin,1) = datasample(noisesample_randsample,timewin);
            end
        end
    case 'randsample_multi'
        noisesample_randsample_multi = noisesample(:);
        if length(noisesample_randsample_multi)<length(signal)
            noisesample_randsample_multi = repmat(noisesample_randsample_multi,ceil(length(signal)/length(noisesample_randsample_multi)),1); %not tested
        end
        timewin = ceil(timewin_perc_of_noisetemplate*Lns);
        samps_val = ceil(Ls/timewin);
        noisesample_new = nan(samps_val*timewin,1);
        for samps = 1:samps_val
            if samps==1
                noisesample_new(1:samps*timewin,1) = datasample(noisesample_randsample_multi,timewin);
            else
                noisesample_new((samps-1)*timewin+1:samps*timewin,1) = datasample(noisesample_randsample_multi,timewin);
            end
        end
        
    case 'trim'
        %Jorges idea goes here
end

%trim to signal length
noisesample_new = noisesample_new(1:Ls);

%adjust noisesignal to yield noisepower in dB by changing amplitude
pwr_sig = rms(signal)^2;
if pwr_sig == 0 %so a signal with zero also works
    pwr_sig = 1;
end
pwr_noise = rms(noisesample_new)^2;
%old version
% what_we_want = (pwr_sig/exp(noisepwr/10));
% x=what_we_want/pwr_noise;
%new
SNR0dB = 10*log10(pwr_sig/pwr_noise);
x = 10^((SNR0dB - noisepwr)/20);

noisesample_new = noisesample_new*x;

%Frequency somehow doesnt work again
Y_ds_p = fft(signal);
Y_ds_p_noise = fft(noisesample_new);
%noisepwr = max(abs(real(Y_ds_p_noise)))/max(abs(real(Y_ds_p)));


Y_new_ds = Y_ds_p + Y_ds_p_noise;%.*noisepwr;

Y_time_new = ifft(Y_new_ds,'symmetric');
%Y_time_new = Y_time_new(1:length(signal));

if p==1
    figure
    hold on
    plot(t,signal,'-b')
    plot(t,Y_time_new,'-r')
end

% %Testing
% %comparing old noise freq spectrum with new noise spectrum
% f_ss_noisesample = samplerate*(0:(Lns/2))/Lns;
% f_ss_sig = samplerate*(0:(Ls/2))/Ls;
% %--------------------------------------
% Y_ds_p_noise_orig = fft(noisesample);
% Y_ss_p_noise_orig = Y_ds_p_noise_orig(1:Lns/2+1);
% %--------------------------------------
% Y_ds_p_noiserep_orig = fft(noisesample_rep(1:Ls));
% Y_ss_p_noiserep_orig = Y_ds_p_noiserep_orig(1:Ls/2+1);
% %--------------------------------------
% Y_ds_p_noiserep_multi = fft(noisesample_rep_multi(1:Ls));
% Y_ss_p_noiserep_multi = Y_ds_p_noiserep_multi(1:Ls/2+1);
% %--------------------------------------
% Y_ds_p_noiserandsample = fft(noisesample_randsample(1:Ls));
% Y_ss_p_noiserandsample = Y_ds_p_noiserandsample(1:Ls/2+1);
% %--------------------------------------
% Y_ds_p_noiserandsample_multi = fft(noisesample_randsample_multi(1:Ls));
% Y_ss_p_noiserandsample_multi = Y_ds_p_noiserandsample_multi(1:Ls/2+1);
% %--------------------------------------
% 
% 
% figure
% hold on
% plot(f_ss_noisesample,abs(Y_ss_p_noise_orig),'-b')  %input noise
% plot(f_ss_sig,abs(Y_ss_p_noiserep_orig),'-g')       %simple repmat
% plot(f_ss_sig,abs(Y_ss_p_noiserep_multi),'-r')      %appending channels
% plot(f_ss_sig,abs(Y_ss_p_noiserandsample),'-k')     %timewindoing datasamples
% plot(f_ss_sig,abs(Y_ss_p_noiserandsample_multi),'-m')%timewindoing datasamples

end