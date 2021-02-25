% -------------------------------------------------------
%
%    run_cvar_bench  - script for running noise study on CVAR-Seg pipeline
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
% This script can be used to rerun the study on the CVAR-Seg pipeline
% small deviations are expected since noise is created randomly from
% clinical noise sample
%
% Example Usage:
%
%    run_cvar_bench
%
% Revision history:



close all
clearvars
clc

addpath(genpath('.'))

load('noisesample.mat'); load('patient_sig.mat');
tic

%=====================================
% Setting up all values
%=====================================

%ground truth CV values
numS1                   = 5;
bclS1                   = 600;      %ms
S1cv_l                  = 660;      %ms
S2list                  = [[500:-50:350]';[340:-10:180]'];   %ms
ERPid                   = length(S2list)-1;         % id of S2list

exp_drop                = 1000; %1000 would be remodeled case
exp_curv                = 0.1;%0.1 would be remodeled case
cv_steady               = 650; %650
erp_ms                  = S2list(ERPid);
rest_curve              = @(exp_drop,exp_curv,s2,erp_cv,cv_steady) -exp_drop*exp(-exp_curv*(s2-erp_cv))+cv_steady;
figure
hold on
plot(S2list,rest_curve(exp_drop,exp_curv,S2list,erp_ms,cv_steady),'-*')
cv_l_list               = rest_curve(exp_drop,exp_curv,S2list,erp_ms,cv_steady); %[[800:-10:660]';[650:-90:200]'];   % conduction velocity in fiber orientation [mm/s]
if sum(cv_l_list<0) > numel(S2list)-ERPid+1
    fprintf('CV should not be negative or 0 before ERP id\n') 
end
anisotropy_ratio        = 1.5;
cv_t_list               = cv_l_list/anisotropy_ratio;        % conduction velocity vertical to fiber orientation [mm/s]

%amplitdue restitution
amp_steady              = 2;
exp_drop_amp            = 2; %1000 would be remodeled case
exp_curv_amp            = 0.05;%0.1 would be remodeled case
atrial_amp_lst          = rest_curve(exp_drop_amp,exp_curv_amp,S2list,erp_ms,amp_steady);
if sum(atrial_amp_lst<=0) > numel(S2list)-ERPid+1
    fprintf('AMP should not be negative or 0 before ERP id\n') 
end
figure
hold on
plot(S2list,atrial_amp_lst,'-*')
ylim([0 5])

%defining default geometrical values of Lasso
def_cat_rad             = 10;         % radius of catheter in mm
numofelecs              = 10;         % not exchangable for now
samplerate              = 5000;      % samplerate at which signal is put together
samplerate_down         = 1000;       % samplerate signals is downsample to
%defining CS
numofelecsCS            = 4; 
dist_step               = 5; %mm
dist_2_stim_elec        = 30:dist_step:30+(numofelecsCS-1)*dist_step;
cv_2_cs_elecs           = cv_l_list; %mm/s

%S1S2 protocol definitions
s2s1dist                = 2*bclS1; %1.4*bclS1;
buffer                  = 1800;      % ms
ERPid;                  %->check CV section

%single atrial activity
s2_aa_sig_width         = 50;      % ms %Not really used - ToDo

%singlestim
t_all                   = 8;
t1_proz                 = 0.25;         %proz of time of first to second part of stimpulse
S                       = 4.7;            %amp in mV
s2                      = 1;            %proz of amplitude of second part compared to first part
tau_edge                = 0.01;
tau_plateau             = 4.5;
B                       = 0;
negbuff                 = t_all*t1_proz;

%noise
dB_noise = [40:-5:10 9:-1:1 0:-0.5:-10];%[40:-5:5 4:-1:1 0:-1:-10];%[80:-5:40 39:-1:10]; %68 for comparison with clin sig
samplerate_noise_seg = 953.674; 

%gaussian noise
noise_constant = 0.005;%0.005; %amp of noise in dB?
baseline_wander = 1;%1; %bool
baseline_mean_amp = 2; %mV
interstitial_beats = 0; %amount of interstitial beats

%settings
save_pics = 0;
ploton = 0;



%====================================================================
% Creating synthetic signal traces
%====================================================================

[LATS1,el_dist,el_pos,~] = getLATtheo(def_cat_rad,S1cv_l,S1cv_l/anisotropy_ratio);
LATS1 = LATS1*1000;  %function gives value in s -> we need ms
for i=1:size(cv_l_list,1)
    LATS2list(:,i) = getLATtheo(def_cat_rad,cv_l_list(i),cv_t_list(i)) * 1000; %function gives value in s
end
LATS1 = LATS1(1:end-1,:); %
LATS2list = LATS2list(1:end-1,:); %
el_dist = el_dist(1:end-1);

for i=1:size(cv_2_cs_elecs,1)
    LATS2list_CS(:,i) = dist_2_stim_elec/cv_2_cs_elecs(i) * 1000;
end
LATS1_CS = LATS2list_CS(:,1);

%=====================
% Creating Templates
%=====================
stim_temp      = createStimTemplate(t1_proz,t_all,S,s2,tau_edge,tau_plateau,B,samplerate); %in ms
stim_temp(:,2) = -stim_temp(:,2); %flipped
f_stim = figure;
save_name = 'f_stim';
plot(stim_temp(:,1),stim_temp(:,2),'-b','LineWidth',2)
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
xlim([-2 6.1])

interstitial_sig_temp = createBipolarAtrialSignals(20,samplerate,'valley');
atrial_sig_temp = createBipolarAtrialSignals(20,samplerate,'downup');
atrial_sig_temp(:,2) = flip(atrial_sig_temp(:,2));


%=====================
% Creating Signal
%=====================
[sigTrace, t_list, logicS2, t_aa_s2, t_aa_s1, sigTraceCS] = createSignalTrace(numS1,bclS1,s2s1dist,...
                LATS1,S2list,LATS2list,stim_temp,atrial_sig_temp,atrial_amp_lst,samplerate,ERPid,buffer,negbuff,...
                LATS1_CS,LATS2list_CS, interstitial_sig_temp, interstitial_beats);
sigTraceCS  = sigTraceCS(1:4); %since CS only has 4 bipolar elecs

ECG = createECGTrace(t_list,sigTrace{1}(:,1));


%===============
% DOWNSAMPLING
%===============
x_newsample = [sigTrace{1}(1,1): 1000/(samplerate_down) :sigTrace{1}(end,1) ]';

x_current = sigTrace{1}(:,1);
mat_trace = cell2mat(sigTrace);
sig_in_sp =  mat_trace(:,2:2:end); clearvars mat_trace;%sigTrace_noise;%
sig_in_cs = sigTraceCS;
for el = 1:size(sig_in_sp,2)
    sigTrace_downsamp(:,el) = interp1(x_current,sig_in_sp(:,el),x_newsample);
end
for el_cs = 1:size(sig_in_cs,2)
    sigTrace_downsamp_CS(:,el_cs) = interp1(sig_in_cs{el_cs}(:,1),sig_in_cs{el_cs}(:,2),x_newsample);
end

ecg_Trace_uniform(:,1) = interp1(ECG(:,1),ECG(:,2),x_newsample);

correl=[];
%%
for noiselvl_i = 1:length(dB_noise)
    %===============
    % ADDING NOISE
    %===============
    
    mat_trace = cell2mat(sigTrace);
    sig_in = mat_trace(:,2:2:end);
    samplerate_down = 1/(mat_trace(2,1)-mat_trace(1,1))*1000;samplerate; %in case change in samplerate is wanted
    clearvars mat_trace;
    sig_in = sigTrace_downsamp;
    sig_in_CS = sigTrace_downsamp_CS;
    samplerate_down = 1/(x_newsample(2,1)-x_newsample(1,1))*1000;
    if noise_constant~=0
        for el = 1:size(sig_in_sp,2)
            sigTrace_noise(:,el)    = addnoise_noisetemp(sig_in(:,el),samplerate_down,Tt,samplerate_noise_seg,dB_noise(noiselvl_i),'randsample_multi',0);
        end
        for el = 1:size(sig_in_CS,2)
            sigTrace_noise_CS(:,el) = addnoise_noisetemp(sig_in_CS(:,el),samplerate_down,Tt,samplerate_noise_seg,dB_noise(noiselvl_i),'randsample_multi',0);
        end
    end
    if baseline_wander ~= 0
        sigTrace_noise    = add_baselinewander(sigTrace_noise,samplerate_down,baseline_mean_amp,1);
        sigTrace_noise_CS = add_baselinewander(sigTrace_noise_CS,samplerate_down,baseline_mean_amp,1);
    end
    
    
    %===============
    % FILTERING
    %===============
    
    sigTrace_final = ECG_High_Low_Filter(sigTrace_noise,samplerate_down,3,300);
    sigTrace_final_CS = ECG_High_Low_Filter(sigTrace_noise_CS,samplerate_down,3,300);
        
    %========================Paper Plot ========
    f_temp = figure;
    save_name = 'f_sig_filt';
    hold on
    stimsig_cut = sigTrace_final(74652:74378+425,6);
    aa_cut = (1.7*sigTrace_final(74378+416:75057,6));
    plot([stimsig_cut; aa_cut],'-r','LineWidth',2);
    plot(5*patient_sig(:),'-b','LineWidth',2);
    xlim([0 400])
    xlabel('samples')
    ylabel('amplitude (a.u)')
    legend('synthetic signal','clinical signal')
    box on
    tcor = corrcoef([stimsig_cut; aa_cut(1:end-15)],patient_sig(:));
    correl(end+1) = tcor(2,1);
  
    %run signal segementation & detection script
    %CS Data is ammended as ones (dummydata)
    ydata = [sigTrace_final_CS sigTrace_final];%[ones(size(sigTrace_final,1),4) sigTrace_final];
    leads = {'CS';'CS';'CS';'CS';'SPI';'SPI';'SPI';'SPI';'SPI';'SPI';'SPI';'SPI';'SPI'};
    elecs = {'1-2';'3-4';'5-6';'7-8';'1-2';'2-3';'3-4';'4-5';'5-6';'6-7';'7-8';'8-9';'9-10'};
    
    xyzin = permute(repmat([rand(4,3) ; [el_pos , zeros(1,size(el_pos,1))'] ],[1 1 25]),[3 2 1]);
    loc_in = location('benchmark',samplerate_down,xyzin,elecs,leads);
    ECGsig = signal('ydata',ecg_Trace_uniform);%for now simple approach
    
    %=========================================
    %run CVAR-seg script: s1s2peakdetect
    %=========================================
    a=s1s2peakdetect(...
        ydata,...
        samplerate_down,...
        leads,...
        elecs,...
        'bipolar','name','BE10_anteriorSPI12_2020_01_28',...
        's1time',600,'s2times',[[500:-50:350] [340:-10:180]],'ERP',180,...
        'map',map(),...
        'location',loc_in,...
        'ECG',ECGsig);
    a.start_detection
    
    %check benchmark against result of script
    result.nums1inblock     = [numS1 a.numof_S1_in_one_block numS1-a.numof_S1_in_one_block];
    result.numofblocks      = [size(S2list,1) a.S1_blocknum size(S2list,1)-a.S1_blocknum];
    result.bcl_ms           = [bclS1 a.S1_time/a.samplerate*1000 bclS1-a.S1_time/a.samplerate*1000];
    result.erp_ms           = [S2list(ERPid) a.s2_steps(a.erp_id)/a.samplerate*1000 S2list(ERPid)-a.s2_steps(a.erp_id)/a.samplerate*1000];
    result.s2aa_sigwidth_ms = [s2_aa_sig_width a.guess_siglength_stim/a.samplerate*1000 s2_aa_sig_width-a.guess_siglength_stim/a.samplerate*1000];
    
    result.s2dt_ms          = {S2list a.s2_steps'/a.samplerate*1000 S2list(1:a.erp_id)-a.s2_steps(1:a.erp_id)'/a.samplerate*1000};
    
    startstopmid            = 3;
    s2p_script              = nan(size(LATS1,1),size(S2list,1));
    tmp                     = a.timeline.getS2list_chan_all;
    s2p_script(:,1:size(tmp,2))= tmp(:,:,startstopmid)./a.samplerate.*1000;
    s2p_orig                = repmat(t_list(logical(logicS2)),1, size(LATS1,1))';
    s2p_orig(:,ERPid:end)   = nan;
    result.s2pos_ms         = [{s2p_orig}   {tmp(:,:,startstopmid)./a.samplerate.*1000}      {s2p_orig - s2p_script}];
    
    s2aaloc_script          = nan(size(t_aa_s2'));
    s2aaloc_script(:,1:size(a.atrial_segments.s2_aa_loc_of_aa,2))= a.atrial_segments.s2_aa_loc_of_aa./a.samplerate*1000;                        
    result.aapos_ms         = [{t_aa_s2'}   {a.atrial_segments.s2_aa_loc_of_aa./a.samplerate*1000}  {t_aa_s2' - s2aaloc_script}]; %ms
    
    s2aadist_ms_script      = nan(size(LATS2list));
    s2aadist_ms_script(:,1:size(a.atrial_segments.s2_aa_dist_fs,2))= a.atrial_segments.s2_aa_dist_fs/a.samplerate*1000;                        
    result.LATs_ms          = [{LATS2list}  {a.atrial_segments.s2_aa_dist_fs/a.samplerate*1000}     {LATS2list - s2aadist_ms_script}];

    s2aap2p_ms_script       = nan(size(LATS2list));
    s2aap2p_ms_script(:,1:size(a.signal_properties.s2_aa_p2p,2))= a.signal_properties.s2_aa_p2p;                        
    theor_amps = abs(max(atrial_sig_temp(:,2))*atrial_amp_lst-min(atrial_sig_temp(:,2))*atrial_amp_lst);
    reform_s2aaamp          = [repmat(theor_amps(1:ERPid-1)' ,size(LATS2list,1),1) nan(size(a.signal_properties.s2_aa_p2p,1),size(S2list,1)-ERPid+1)];%[repmat(s2_aa_amp,size(a.signal_properties.s2_aa_p2p,1),ERPid-1)  zeros(size(a.signal_properties.s2_aa_p2p,1),size(S2list,1)-ERPid+1)];
    result.amps2aa_mV       = [{reform_s2aaamp} {a.signal_properties.s2_aa_p2p} {reform_s2aaamp - s2aap2p_ms_script}]; %mV
    
    s2aainvLAT_ms_script    = nan(size(LATS2list));
    s2aainvLAT_ms_script(:,1:size(a.atrial_segments.s2_aa_dist_fs,2))= a.atrial_segments.s2_aa_dist_fs/a.samplerate*1000;                        
    result.invLATs_ms       = [{1 ./ LATS2list}  {1 ./ (a.atrial_segments.s2_aa_dist_fs/a.samplerate*1000)}     {((1 ./ LATS2list - 1 ./ (s2aainvLAT_ms_script)))}];
    
    s2aaCV_ms_script        = nan(size(LATS2list));
    s2aaCV_ms_script(:,1:size(result.invLATs_ms{1,2},2))= result.invLATs_ms{1,2}/1000;                        
    result.cv_mmps          = [{el_dist./(LATS2list/1000)} {el_dist./(result.invLATs_ms{1,2}/1000)}  {el_dist.*result.invLATs_ms{1,3}} ]; % {el_dist./(LATS2list/1000) - el_dist./(s2aaCV_ms_script) }
    
    overall_result(noiselvl_i) = result;
    overall_sig(noiselvl_i,:,:) = ydata;
    
    fn = fieldnames(result);
    er = [];
    for fi = 1:size(fn,1)
        if iscell(result.(fn{fi}))
            err.noise = dB_noise(noiselvl_i);
            
            %Median Error
            for chani = 1:size(result.(fn{fi}){1,3},1)
                %mean/med of single channel over all s2
                err.(fn{fi}).chan.mean   = mean(result.(fn{fi}){1,3},2,'omitnan');
                err.(fn{fi}).chan.median = median(result.(fn{fi}){1,3},2,'omitnan');
                %mean/med of single s2 over all channels
                err.(fn{fi}).s2.mean   = mean(result.(fn{fi}){1,3},1,'omitnan');
                err.(fn{fi}).s2.median = median(result.(fn{fi}){1,3},1,'omitnan');
            end
            err.(fn{fi}).meas.mean	 = mean(result.(fn{fi}){1,3},'omitnan');
            err.(fn{fi}).meas.median = median(result.(fn{fi}){1,3},'omitnan');
        end
    end
    overall_error(noiselvl_i) = err; %errors for a single noise lvl
    close all
    
    if ploton
        for fi = 1:size(fn,1)
            if iscell(result.(fn{fi}))
                f = figure;
                set(f,'Position',[100 100 1200 600])
                heatmap(round(result.(fn{fi}){1,3},1))%heatmap(abs(result.(fn{fi}){1,3}))
                unitstr = fn{fi}(end-1:end);
                stri = ['Error of ' fn{fi} ' compared to synthetic signal ' unitstr];
                title(stri)
                xlabel('S2 beat')
                ylabel('Electrodes')
                colorbar
                save_name = [fn{fi} ' error syntheticsig samplerate ' num2str(samplerate_down) ' noise ' num2str(noise_constant) ' baselw' num2str(baseline_wander)];
                if save_pics==1
                    savefig(f,[save_name '.fig'])
                    print(f,[save_name '.png'],'-dpng','-r300')
                end
                
            elseif size(result.(fn{fi}),3) == 1
                fprintf('Error of %s : %f \n',fn{fi},result.(fn{fi})(1,3))
            else %is never used for now
                f = figure;
                set(f,'Position',[100 100 1200 600])
                unitstr = fn{fi}(end-1:end);
                stri = ['Error of ' fn{fi} ' compared to synthetic signal ' unitstr];
                plot(1:size(result.(fn{fi}),1),abs(result.(fn{fi})(:,3)),'-xb','LineWidth',2)
                title(stri)
                xlabel('Electrodes')
                ylabel('Error (ms)')
                %xlim([1 size(result.(fn{fi}),1)])
                save_name = [fn{fi} ' error syntheticsig samplerate ' num2str(samplerate_down) ' noise ' num2str(noise_constant) ' baselw' num2str(baseline_wander)];
                if save_pics==1
                    savefig(f,[save_name '.fig'])
                    print(f,[save_name '.png'],'-dpng','-r300')
                end
            end
        end
        
    end
    
    
end

plot_results
