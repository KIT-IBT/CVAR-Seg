% -------------------------------------------------------
%
%    filter_stimartifact_matchedfilt - matched filter to reduce stimulationa artifact 
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

function [filt_signal] = filter_stimartifact_matchedfilt(signal,stimtemplate,upsample,scale2peak)
%FILTER_STIMARTIFACT_MATCHEDFILT:
%Uses the matched filter approach to subtract stimtemplate from the signal
%upsample is the factor by wich the signal is upsampled
%scale2peak is integer and describes wich peak the signal should be scaled
%to -> If a stimulus as several decreasing amplitudeslopes 2, would scale
%to the second peak in time.

%TODO check input
filt_signal = [];

stp_up = interp(stimtemplate,upsample);
s_up = interp(signal,upsample);

sz_dif = abs(length(stp_up) - length(s_up));
stp_up = [stp_up;zeros(sz_dif,1)];
% if mod(sz_dif,2)==0
%     stp_up = [zeros(sz_dif/2,1); stp_up; zeros(sz_dif/2,1)];
% else
%     stp_up = [ceil(zeros(sz_dif/2,1)); stp_up; floor(zeros(sz_dif/2,1))];
% end

diff_sig = abs(diff(s_up)); %worked well without abs()
%[~,diff_loc]=max(diff_sig);
[s_up_diff_val, s_up_diff_loc] = findpeaks( diff_sig );
pks_sort_diff = sortrows([s_up_diff_val s_up_diff_loc],'descend');
buffer=20;
peakcoice_diff = 1; %1 with buffer is better than 2 without buffer!
if pks_sort_diff(peakcoice_diff,2)+buffer > length(s_up)
    s_up_diff_cut = s_up(1:pks_sort_diff(peakcoice_diff,2));
else
    s_up_diff_cut = s_up(1:pks_sort_diff(peakcoice_diff,2)+buffer);
end
[val , loc]=xcorr(s_up_diff_cut,stp_up);
[m_val , m_loc]=max(val);

%[val , loc]=xcorr(s_up,stp_up);
%[m_val , m_loc]=max(val);
if m_loc < 1
    %TODO (just switch sig with temp?)
    fprintf('TODO!\n')
end

%find signal amplitude of wanted peak
[s_val, s_loc] = findpeaks( s_up.^2 );
pks_sort_s = sortrows([sqrt(s_val) s_loc],'descend');
[~ , pk_loc_scale]= max(pks_sort_s(1:scale2peak,2));%from first 2 peaks with highest Amp chose the one closer to stimulus so it will be deleted better
pk_choice_s = [s_up(pks_sort_s(pk_loc_scale,2)) pks_sort_s(pk_loc_scale,2)];
s_amp = pk_choice_s(1);


%find template amplitude of wanted peak
[stp_val, stp_loc] = findpeaks( stp_up.^2 );
pks_sort_stp = sortrows([sqrt(stp_val) stp_loc],'descend');
[~ , pk_loc_scale]= max(pks_sort_stp(1:scale2peak,2));%from first 2 peaks with highest Amp chose the one closer to stimulus so it will be deleted better
pk_choice_stp = [stp_up(pks_sort_stp(pk_loc_scale,2)) pks_sort_stp(pk_loc_scale,2)];
stp_amp = pk_choice_stp(1);

scaleingfactor = s_amp/stp_amp;
if scaleingfactor<0
    scaleingfactor=1;
end
%subtract
stp_shift_up = circshift(stp_up,m_loc)*(scaleingfactor);

% if m_loc>1
%     stp_up_shift = [zeros(m_loc-length(s_up),1); stp_up];
% else
%     stp_up_shift = [stp_up;zeros(m_loc-length(s_up),1)];
% end
% stp_up_shift = stp_up_shift(1:length(s_up))*(scaleingfactor);
s_sub = s_up - stp_shift_up;


%downsample
s_sub_down = decimate(s_sub,upsample);
stp_shift_down = decimate(stp_shift_up,upsample);

filt_signal = s_sub_down;

%===========
%check plots
%===========
%downsampled
% figure
% hold on
% plot(signal,'b-')
% plot(stp_shift_down,'b-')
% 
% figure
% ax1 = subplot(3,1,1);
% plot(signal,'b-')
% ax2 = subplot(3,1,2);
% plot(stp_shift_down,'b-')
% ax3 = subplot(3,1,3);
% plot(s_sub_down,'b-')

% %upsampled
% figure
% hold on
% plot(s_up,'b-')
% plot(stp_shift_up,'b-')

% figure
% ax1 = subplot(3,1,1);
% plot(s_up,'b-')
% ax2 = subplot(3,1,2);
% plot(stp_shift_up,'r-')
% ax3 = subplot(3,1,3);
% plot(s_sub,'g-')
% linkaxes([ax1 ax2 ax3],'x');

end

