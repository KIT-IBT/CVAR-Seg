%%
%Amps after downsampling
amp_timewin = 10; %samples
amp_thresh = max(theor_amps)/2;
LATS2AAlist = [nan(size(t_aa_s2,1),4), t_aa_s2];%(LATS2list*samplerate_down)+buffer+negbuff;%
LATS2AAlist(LATS2AAlist==0)=nan;
probs = [0.25 0.5 .75];
aamp_downsampled_err = nan(length(dB_noise),size(overall_sig,3),size(LATS2AAlist,2));
for noiselvl_i = 1:size(overall_sig,1)
    for chan = 1:size(overall_sig,3)
        for ilats = 1:size(LATS2AAlist,1)
            if isnan(LATS2AAlist(ilats,chan))
                continue
            else
                il_start = round(LATS2AAlist(ilats,chan) -amp_timewin);
                il_end = round(LATS2AAlist(ilats,chan) +amp_timewin*2);
                try
                    [aamp_downsampled_min(noiselvl_i,ilats,chan),loc_min_amp] = min(overall_sig(noiselvl_i,round(LATS2AAlist(ilats,chan))+2:il_end,chan));
                catch
                    aamp_downsampled_min(noiselvl_i,ilats,chan) = nan;
                end
                try
                    [pks,locs] = findpeaks(overall_sig(noiselvl_i,round(LATS2AAlist(ilats,chan))+loc_min_amp-1-amp_timewin:round(LATS2AAlist(ilats,chan))+loc_min_amp-1,chan));
                    pks_sort = sortrows([pks' locs'],'descend');
                    aamp_downsampled_max(noiselvl_i,ilats,chan) = pks_sort(end,1);
                catch
                    aamp_downsampled_max(noiselvl_i,ilats,chan) = nan;
                end
                
                aamp_downsampled(noiselvl_i,ilats,chan) = aamp_downsampled_max(noiselvl_i,ilats,chan)-aamp_downsampled_min(noiselvl_i,ilats,chan);
            end
            
            %new err calc
            aamp_downsampled_err(noiselvl_i,ilats,chan) = aamp_downsampled(noiselvl_i,ilats,chan) - overall_result(noiselvl_i).amps2aa_mV{2}(chan-4,ilats);
        end
    end
    m_amp_er(noiselvl_i) = median(aamp_downsampled_err(noiselvl_i,:,:),'all');
    quants_amp_new(noiselvl_i,:) = quantile(aamp_downsampled_err(noiselvl_i,:,:),probs,'all')';
end

%%
exclude_stimchan = 1; % <-- either 0 or chanid of stimchan
allids = 1:size(overall_result(1).aapos_ms{1,3},1);
evalids = allids(~ismember(allids,exclude_stimchan));
for ii=1:size(overall_error,2)
    if exclude_stimchan
        ov_aapos_med_err(ii) = median(overall_result(ii).aapos_ms{1,3}(evalids,:),'all','omitnan');
        ov_aaamp_med_err(ii) = median(overall_result(ii).amps2aa_mV{1,3}(evalids,:),'all','omitnan');
    else
        ov_aapos_med_err(ii) = median(overall_result(ii).aapos_ms{1,3}(:),'omitnan');
        ov_aaamp_med_err(ii) = median(overall_result(ii).amps2aa_mV{1,3}(:),'omitnan');
        %ov_aaamp_mean_err(ii) = mean(overall_result(ii).amps2aa_mV{1,3}(:),'omitnan');
        ov_meas_max_err(ii) = max(overall_error(ii).aapos_ms.meas.median);
    end
end
%%
%Quantiles % missdetect counts
fnq = fieldnames(overall_error);
mean_error_2_high_SNR_case = quantile((overall_result(1).(fnq{4}){1,3}(end-1,:)),probs,'all')';
for ii=1:size(overall_error,2) %dB
    for filds = 2:size(fnq,1)
        probs = [0.25 0.5 .75];
        if exclude_stimchan
            quants(ii,filds,:) = quantile((overall_result(ii).(fnq{filds}){1,3}(evalids,:)),probs,'all')'; %abs verhindert das 75 quantil größer als 25 quantil
            missdect(ii,filds,:) = sum(isnan(overall_result(ii).(fnq{filds}){1,3}(evalids,:)),'all') ...
                - (size(overall_result(ii).(fnq{filds}){1,3},2)-ERPid+1)*size(overall_result(ii).(fnq{filds}){1,3},1);
        else
            quants(ii,filds,:) = quantile((overall_result(ii).(fnq{filds}){1,3}(:)),probs); %abs verhindert das 75 quantil größer als 25 quantil
            missdect(ii,filds,:) = sum(isnan(overall_result(ii).(fnq{filds}){1,3}(:))) ...
                - (size(overall_result(ii).(fnq{filds}){1,3},2)-ERPid+1)*size(overall_result(ii).(fnq{filds}){1,3},1);
        end
    end
end
%%
%S1 time detection

% %AMP ERROR OVER DB ==LEGACY==
% f_medAAAmperror = figure;
% hold on
% plot(dB_noise(1:size(ov_aaamp_med_err,2)),quants(:,6,1),'-b')
% plot(dB_noise(1:size(ov_aaamp_med_err,2)),quants(:,6,2),'-k')
% plot(dB_noise(1:size(ov_aaamp_med_err,2)),quants(:,6,3),'-r')


f_medAAAmperror = figure;
hold on
plot(dB_noise(1:size(ov_aaamp_med_err,2))',quants_amp_new(:,2),'-ok','LineWidth',1.0)
shadedErrorBar(dB_noise(1:size(ov_aaamp_med_err,2))',quants_amp_new(:,2),[abs(quants_amp_new(:,2)-quants_amp_new(:,3))';zeros(size(quants_amp_new(:,2)))'],'lineProps','-r')%upperbound
shadedErrorBar(dB_noise(1:size(ov_aaamp_med_err,2))',quants_amp_new(:,2),[zeros(size(quants_amp_new(:,2)))';abs(quants_amp_new(:,2)-quants_amp_new(:,1))'],'lineProps','-b')
ylabel('Amplitude error (mV)')
xlabel('SNR (dB)')
legend({['Q2'],['Q3'],['Q1']})
ylim([-0.5 0.5])
set(gca,'FontSize',12)
yticks(-0.5:0.1:0.5)
box on
grid on
save_name = 'AA Amp error quartiles';
% savefig(f_medAAAmperror,['Pics/' save_name '.fig'])
% print(f_medAAAmperror,['Pics/' save_name],'-dpng','-r300')
% hold on

aamp_perc = aamp_downsampled_err./aamp_downsampled*100;
quants_amp_new_perc = quantile(aamp_perc,probs,[2 3]);
f_medAAAmperror_perc = figure;
hold on
plot(dB_noise(1:size(ov_aaamp_med_err,2))',quants_amp_new_perc(:,2),'-ok','LineWidth',1.5)
shadedErrorBar(dB_noise(1:size(ov_aaamp_med_err,2))',quants_amp_new_perc(:,2),[abs(quants_amp_new_perc(:,2)-quants_amp_new_perc(:,3))';zeros(size(quants_amp_new_perc(:,2)))'],'lineProps','-r')%upperbound
shadedErrorBar(dB_noise(1:size(ov_aaamp_med_err,2))',quants_amp_new_perc(:,2),[zeros(size(quants_amp_new_perc(:,2)))';abs(quants_amp_new_perc(:,2)-quants_amp_new_perc(:,1))'],'lineProps','-b')
ylabel('Amplitude error (%)')
xlabel('SNR (dB)')
legend({['Q2'],['Q3'],['Q1']})
set(gca,'FontSize',12)
box on
grid on
save_name = 'AA Amp error percent';
% savefig(f_medAAAmperror_perc,['Pics/' save_name '.fig'])
% print(f_medAAAmperror_perc,['Pics/' save_name],'-dpng','-r300')
% hold on



brAMP_col=[];brAMP_row=[];
for noiselvl_i = 1:size(overall_result,2)
    tmp = aamp_downsampled_err(noiselvl_i,:,:);
    tmp(tmp==0)= nan;
    brAMP(:,noiselvl_i) = tmp(:);
    brAMP_col = [brAMP_col;tmp(:)];
    brAMP_row = [brAMP_row;ones(length(tmp(:)),1)*noiselvl_i];
end


%LAT ERROR OOVER DB
f_medLATerror = figure;
hold on
plot(dB_noise(1:size(ov_aapos_med_err,2))',quants(:,4,2),'-ok', 'LineWidth',1.0)
shadedErrorBar(dB_noise(1:size(ov_aapos_med_err,2))',quants(:,4,2),[abs(quants(:,4,2)-quants(:,4,3))';zeros(size(quants(:,4,2)))'],'lineProps','-r')%upperbound
shadedErrorBar(dB_noise(1:size(ov_aapos_med_err,2))',quants(:,4,2),[zeros(size(quants(:,4,2)))';abs(quants(:,4,2)-quants(:,4,1))'],'lineProps','-b')
xlabel('SNR (dB)')
legend({['Q2'],['Q3'],['Q1']})
ylabel('LAT error (ms)')
xlabel('SNR (dB)')
ylim([-3 3])
xlim([-7 40])
set(gca,'FontSize',12)
box on
grid on
save_name = 'AA LAT error quartiles';
% savefig(f_medLATerror,['Pics/' save_name '.fig'])
% print(f_medLATerror,['Pics/' save_name],'-dpng','-r300')

brLAT_col=[];brLAT_row=[];
for noiselvl_i = 1:size(overall_result,2)
    tmp = overall_result(noiselvl_i).(fnq{4}){1,3}(evalids,:);
    brLAT(:,noiselvl_i) = tmp(:);
    brLAT_col = [brLAT_col;tmp(:)];
    brLAT_row = [brLAT_row;ones(length(tmp(:)),1)*noiselvl_i];
end

%Num of NAN values of LAT(missdect(:,6,:)) & AMP(missdect(:,4,:))
f_mLAT = figure;
plot(dB_noise(1:size(ov_aapos_med_err,2))', missdect(:,4,:)', '-*k', 'LineWidth',2)
hold on
%title('Number of not-detected LATs')
ylabel('num of missed LATs (#)')
xlabel('SNR (dB)')
box on
grid on
set(gca,'FontSize',12)
save_name = 'missedLATs';
% savefig(f_mLAT,['Pics/' save_name '.fig'])
% print(f_mLAT,['Pics/' save_name],'-dpng','-r300')

f_mLAT_perc = figure;
plot(dB_noise(1:size(ov_aapos_med_err,2))', missdect(:,4,:)/((size(overall_result(ii).(fnq{filds}){1,3},1)*(ERPid-1))-(ERPid-1)), '-ok', 'LineWidth',1.5)
hold on
%title('Number of not-detected LATs')
ylabel('missed LATs (%)')
xlabel('SNR (dB)')
set(gca,'FontSize',12)
xlim([-7 40])
box on
grid on
save_name = 'missedLATsperc';
% savefig(f_mLAT_perc,['Pics/' save_name '.fig'])
% print(f_mLAT_perc,['Pics/' save_name],'-dpng','-r300')


%heatmaps
f_hLAT = figure;
hidx = 2;
hstr = 'aapos_ms';'amps2aa_mV';
set(f_hLAT,'Position',[100 100 1200 600])
heatmap(round(overall_result(hidx).(hstr){1,3},1))%heatmap(abs(result.(fn{fi}){1,3}))
unitstr = fn{fi}(end-1:end);
stri = ['Error of atrial LAT compared to synthetic signal (ms)'];
title(stri)
xlabel('S2 beat (#)')
ylabel('Bipolar Channels')
colorbar
% save_name = [hstr ' error syntheticsig samplerate ' num2str(samplerate_down) ' noise ' num2str(dB_noise(noiselvl_i)) ' baselw' num2str(baseline_wander)];
% savefig(f_hLAT,['Pics/' save_name '.fig'])
% print(f_hLAT,['Pics/' save_name],'-dpng','-r300')


%plot simple cv restitution curves (single measurement)
%%
cv_err=[];
hstr = 'cv_mmps';
rf = figure;
hold on; box on; grid on;
l={};
%plot(S2list(1:ERPid-1),cv_l_list(1:ERPid-1),'-*','LineWidth',2)
for hidx = 1%:10%length(dB_noise)
    for r = [3 4 5 6 7 8 9]%1:size(overall_result(hidx).(hstr){1, 2},1)
        plot(S2list(1:ERPid-1),overall_result(hidx).(hstr){1,1}(r,1:ERPid-1),'-*b','LineWidth',2)
        l{end+1} = num2str(r);
        erp_cv_id = find(S2list==overall_result(hidx).erp_ms(1,2))-1;
        if erp_cv_id>size(overall_result(hidx).invLATs_ms{1,2},2);erp_cv_id=size(overall_result(hidx).invLATs_ms{1,2},2);end
        try
            x = el_dist(r)./(overall_result(hidx).LATs_ms{1,2}(r,1:erp_cv_id))*1000;
            inn = ~isnan(x);
            cvintp0 = interp1(S2list(inn),x(inn),'linear','pp');
            cvintp = fnval(cvintp0,S2list(1:erp_cv_id));
        catch
            cvintp = [];
        end
        if numel(cvintp)>1
            plot(S2list(1:erp_cv_id), cvintp,'-*r','LineWidth',1.5)
            cv_err(:,end+1)=overall_result(hidx).(hstr){1,1}(r,1:ERPid-1)'-cvintp;
        end
        l{end+1} = num2str(r);
    end
end
legend(l)
set(gca,'FontSize',12)
ylim([0 800])
xlim([150 max(S2list)])
cv_err_mean = mean(cv_err);
LAT_err_mean = mean(el_dist(5)./(cv_err(:,3))');
% save_name = 'f_CV_rest_nonoise';
% savefig(rf,['Pics/' save_name '.fig'])
% print(rf,['Pics/' save_name],'-dpng','-r300')
%%

%plot simple amp restitution curves (single measurement)
%%
hstr = 'amps2aa_mV';
m2mfac = max(atrial_sig_temp(:,2))/max(abs(atrial_sig_temp(:,2))) + 1; %since we only scale to max peak
rf = figure;
hold on; box on; grid on;
theor_amps = abs(max(atrial_sig_temp(:,2))*atrial_amp_lst-min(atrial_sig_temp(:,2))*atrial_amp_lst);
for hidx = 1%:10%length(dB_noise)
    for chani = 3:size(overall_result(hidx).(hstr){1, 2},1)
        erp_amp_id = find(S2list==overall_result(hidx).erp_ms(1,2))-1;
        if erp_amp_id>size(overall_result(hidx).(hstr){1,2},2);erp_amp_id=size(overall_result(hidx).(hstr){1,2},2);end
        try
            x = overall_result(hidx).(hstr){1,2}(chani,1:erp_amp_id);
            inn = ~isnan(x);
            ampintp0 = interp1(S2list(inn),x(inn),'linear','pp');
            ampintp = fnval(ampintp0,S2list(1:erp_amp_id));
        catch
            ampintp = [];
        end
        if numel(ampintp)>1
            plot(S2list(1:erp_amp_id), ampintp,'-*r','LineWidth',1.5)
        end
        
    end
end
plot(S2list(1:ERPid-1),theor_amps(1:ERPid-1),'-*b','LineWidth',2)
set(gca,'FontSize',12)
% save_name = 'f_amp_rest_nonoise';
% savefig(rf,['Pics/' save_name '.fig'])
% print(rf,['Pics/' save_name],'-dpng','-r300')

%%

fprintf('Done\n')
toc

%plot all noiselevls
f_nois = figure;
chanp = 11;
for noiselvl_i = 1:length(dB_noise)
    ax(noiselvl_i) = subplot(length(dB_noise),1,noiselvl_i);
    plot((0:size(overall_sig,2)-1)/samplerate_down,squeeze(overall_sig(noiselvl_i,:,chanp)))
    box on
    title([num2str(dB_noise(noiselvl_i)) ' (dB)'] )
    ylabel('Amplitude (a.u.)' )
end
xlabel('Time (s)')
linkaxes(ax,'xy');
xlim([1.5 5.0])
ylim([-5 5])

%Noiselevels Paper plot
n=[1 7 17 25 33];
f_nois = figure;
chanp = 11;
for noiselvl_i = 1:length(n)
ax(noiselvl_i) = subplot(length(n),1,noiselvl_i);
plot((0:size(overall_sig,2)-1)/samplerate_down,squeeze(overall_sig(n(noiselvl_i),:,chanp)))
box on
title([num2str(dB_noise(n(noiselvl_i))) ' (dB)'] )
if noiselvl_i<length(n)
set(gca,'xtick',[])
set(gca,'xticklabel',[])
end
end
linkaxes(ax,'xy');

% save_name = 'noise_lvls';
% savefig(f_nois,['Pics/' save_name '.fig'])
% print(f_nois,['Pics/' save_name],'-dpng','-r300')
            
%%
%save('noise_study_paper_awesome_V4.mat')