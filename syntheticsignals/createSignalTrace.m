% -------------------------------------------------------
%
%    createStimTemplate - puts together signal trace of S1S2 protocol based
%    on LATs and templates for CS and roving catheter
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

function [signal,t_start,logicS2,t_aa_s2,t_aa_s1,signal_cs] = createSignalTrace(numS1,bclS1,s2s1dist,LATS1,S2list,LATS2,sig_stim_temp,atrial_sig_temp,atrial_amp_lst,samplerate,ERPid,buffer,neg_buff,LATS1_CS,LATS2_CS,interbeat_sig_temp,interbeats)
%creates startingpointlist tstart in ms. neg_buf is to start away from zero
%in cas normal buffer is 0 (not recommended)

%factor s2 aa 2 s1 aa amplitude - s2 amp mostly a bit lower
%s2aa_fac = 0.8;

t_start = neg_buff + buffer; %not 0 since we need room to center stimulus on first beat, t_start in 
logicS2 = [];
for i=1:size(S2list,1)
    %cnt=0;
    for s1cnt = 1:numS1
        %cnt=cnt+1;
        if s1cnt ~= 1
            t_start(end+1,1)= t_start(end,1) + bclS1;
        end
        logicS2(end+1,1) = false;
    end
    
    %s2 beat
    t_start(end+1,1)= t_start(end,1) + S2list(i);
    logicS2(end+1,1) = true;
    if S2list(i) ~= S2list(end) %last S2 should nor be followed by another s1s2 protocol
        t_start(end+1,1)= t_start(end,1) + s2s1dist; %so startingpoint of next iteration is correct
    end
end


%Set template on startingpoints
% x & y is for roving catheter & x_cs y_cs is CS
cs_factor_stim = 0.025; %factor for cs stim amp compared to atrial
cs_factor_atrial = 1; 
cs_window_delay = 8 / 1000 *samplerate; %delay between atrial s2 and cs s2 

t_aa_s2 = nan(2,2);
t_aa_s1 = nan(2,2);
l_id = find(logicS2==1);
stim_sz = size(sig_stim_temp,1);
atrial_sz = size(atrial_sig_temp,1);
[~,atrial_temp_center_id] = min(abs(atrial_sig_temp(:,1)));
[~,stim_temp_center_id] = min(abs(sig_stim_temp(:,1)));
atrial_temp_center_id = abs(atrial_sig_temp(1,1));

t_end = (t_start(end) + max(LATS1(:)) + s2s1dist + cs_window_delay )/1000;
x = [0:1/samplerate:t_end]'*1000; %is in ms
samps = length(x);

for el = 1:size(LATS1,1)
    cnts2=0;
    y_sp = zeros(samps,1);
    y_cs = zeros(samps,1);
    %x=[];
    %y=[];
    %x_cs = [];
    %y_cs = [];
    for t = 1:size(t_start,1)
        
        %set stimulus template on the times
        interptemp2sig = interp1(t_start(t) + sig_stim_temp(:,1) , sig_stim_temp(:,2) , x);
        y_sp = sum([y_sp interptemp2sig],2,'omitnan');
        
        interptemp2sig = interp1(t_start(t) + sig_stim_temp(:,1), sig_stim_temp(:,2) , x);
        y_cs = sum([y_cs interptemp2sig*cs_factor_stim],2,'omitnan');
        
        if ismember(t,l_id) %add atrial answer fitting to the s2 cv
            cnts2 = cnts2 + 1;
            if cnts2 < ERPid
                lat_dif = round((abs(atrial_temp_center_id) - LATS2(el,cnts2))/1000*samplerate);
                if lat_dif > 0
                    as = [atrial_sig_temp(:,1) [zeros(lat_dif,1); atrial_sig_temp(lat_dif+1:end,2)]];
                else
                    as = atrial_sig_temp;
                end
                interptemp2sig = interp1(t_start(t) + LATS2(el,cnts2) + as(:,1),as(:,2),x);
                y_sp = sum([y_sp atrial_amp_lst(cnts2)*interptemp2sig],2,'omitnan');
                
%                 interptemp2sig = interp1(t_start(t) + LATS2(el,cnts2) + cs_window_delay + atrial_sig_temp(:,1), ...
%                     atrial_sig_temp(:,2)*cs_factor_atrial , x);
                if el <= size(LATS2_CS,1)
                    interptemp2sig = interp1(t_start(t) + LATS2_CS(el,cnts2) + atrial_sig_temp(:,1), atrial_sig_temp(:,2) , x);
                    y_cs = sum([y_cs interptemp2sig*cs_factor_atrial*atrial_amp_lst(cnts2)],2,'omitnan');
                end
                t_aa_s2(cnts2,el) = [t_start(t) + LATS2(el,cnts2)];
            else
                t_aa_s2(cnts2,el) = nan;
            end
        else %add atrail answer fitting to s1 cv
            lat_dif = round((abs(atrial_temp_center_id) - LATS1(el))/1000*samplerate);
            if lat_dif > 0
                as = [atrial_sig_temp(:,1) [zeros(lat_dif,1); atrial_sig_temp(lat_dif+1:end,2)]];
            else
                as = atrial_sig_temp;
            end
            interptemp2sig = interp1(t_start(t) + LATS1(el) + as(:,1),as(:,2),x);
            y_sp = sum([y_sp atrial_amp_lst(1)*interptemp2sig],2,'omitnan'); %functionvalue
            
            t_aa_s1(t-cnts2,el) = [t_start(t) + LATS1(el)];
            
            %interptemp2sig = interp1(t_start(t) + LATS1(el) + cs_window_delay + atrial_sig_temp(:,1), ...
            %                                             atrial_sig_temp(:,2)*cs_factor_atrial , x);
            if el <= size(LATS2_CS,1)
                interptemp2sig = interp1(t_start(t) + LATS1_CS(el) + atrial_sig_temp(:,1), ...
                                    atrial_sig_temp(:,2)*cs_factor_atrial , x);
                y_cs = sum([y_cs interptemp2sig],2,'omitnan');
            end
            y_cs = sum([y_cs interptemp2sig],2,'omitnan');
            
            %x_cs = [x_cs; t_start(t) + LATS1(el) + cs_window_delay + atrial_sig_temp(:,1)];
            %y_cs = [y_cs; atrial_sig_temp(:,2)*3*cs_factor]; 
        end
    end
    
%     %add buffer to start & end of signal
%     x_buf(:,1) = linspace(0,buffer/1000,buffer/1000 * samplerate);
%     x = [x_buf ; x ; x(end) + 1/samplerate + x_buf];
%     y_sp = [repmat(0,[buffer/1000*samplerate 1]) ; y_sp ; repmat(0,[buffer/1000*samplerate 1]) ]; %buffervalues
    signal{el} = [x y_sp];
    
    %x_cs = [x_buf ; x_cs ; x(end) + 1/samplerate + x_buf];
    %y_cs = [repmat(0,[buffer/1000*samplerate 1]) ; y_cs ; repmat(0,[buffer/1000*samplerate 1]) ];
    signal_cs{el} = [x y_cs];
end


%adding interstitial beats to lasso
same_in_all_chan = 1;
intb_amp = 2;
if interbeats ~= 0
    %find timeslots larger than
    t_win_intb = 450; %ms
    fact = 150/t_win_intb; % 150 is smallest ERP / 
    for el = 1:size(LATS1,1)
        timeline(:,el)      = sort([t_start;t_aa_s1(:,el);t_aa_s2(:,el)]);
        timeline_diff(:,el) = diff(timeline(:,el));
        log_win(:,el)       = timeline_diff(:,el)>t_win_intb;
    end
    
    
    if same_in_all_chan
        [fidx] = find(log_win(:,1)==1);
        chosen_t = randsample(timeline(fidx,1),interbeats);
        %choose rand timeinterval (larger than 20 percen so there is room to
        %iniial signal
        interbeat_buff = fact*t_win_intb;
        
        rand_t = interbeat_buff + (t_win_intb-interbeat_buff).*rand(interbeats,1);
        t_interb = chosen_t + rand_t;
        
        %add signal to all channels to the timeline
        for el = 1:size(LATS1,1)
            for tb = 1:size(t_interb,1)
                interptemp2sig = interp1(t_interb(tb) + interbeat_sig_temp(:,1) , intb_amp*interbeat_sig_temp(:,2) , x);
                signal{el}(:,2) = sum([signal{el}(:,2) interptemp2sig],2,'omitnan');
            end
        end
    else
        fprintf('TODO\n')
%         %choose rand starting timeslot
%         [rx,cx] = find(log_win==1);
%         chosen_t = randsample(timeline(fidx),interbeats);
%         %choose rand timeinterval (larger than 20 percen so there is room to
%         %iniial signal
%         interbeat_buff = 0.2*t_win_intb;
%         
%         rand_t = interbeat_buff + (t_win_intb-interbeat_buff).*rand(interbeats,1);
%         t_interb = chosen_t + rand_t;
%         
%         %add signal to all channels to the timeline
%         for tb = 1:size(t_interb,1)
%             interptemp2sig = interp1(t_interb(tb) + interbeat_sig_temp(:,1) , interbeat_sig_temp(:,2) , x);
%             y_sp = sum([y_sp interptemp2sig],2,'omitnan');
%         end
    end
    %figure
    %plot(timeline(1,:),ones(length(timeline(1,:)),1),'*');
    %plot(timeline_diff(1,:),ones(length(timeline_diff(1,:)),1),'*')
    %plot(signal{el}(:,1),signal{el}(:,2))
end



%for el = 1:size(LATS1,1)
%figure
%plot(signal{el}(:,1),signal{el}(:,2))
%end


fprintf('Finished creating signals\n')
end





%CODE SAVE


% t_start = neg_buff + buffer; %not 0 since we need room to center stimulus on first beat, t_start in 
% logicS2 = [];
% for i=1:size(S2list,1)
%     %cnt=0;
%     for s1cnt = 1:numS1
%         %cnt=cnt+1;
%         if s1cnt ~= 1
%             t_start(end+1,1)= t_start(end,1) + bclS1;
%         end
%         logicS2(end+1,1) = false;
%     end
%     
%     %s2 beat
%     t_start(end+1,1)= t_start(end,1) + S2list(i);
%     logicS2(end+1,1) = true;
%     if S2list(i) ~= S2list(end) %last S2 should nor be followed by another s1s2 protocol
%         t_start(end+1,1)= t_start(end,1) + s2s1dist; %so startingpoint of next iteration is correct
%     end
% end
% 
% 
% %Set template on startingpoints
% % x & y is for roving catheter & x_cs y_cs is CS
% cs_factor = 10;
% cs_window_delay = 4 / 1000 *samplerate;
% 
% t_aa_s2 = nan(2,2);
% t_aa_s1 = nan(2,2);
% l_id = find(logicS2==1);
% stim_sz = size(sig_stim_temp,1);
% atrial_sz = size(atrial_sig_temp,1);
% for el = 1:size(LATS1,1)
%     cnts2=0;
%     x=[];
%     y=[];
%     x_cs = [];
%     y_cs = [];
%     for t = 1:size(t_start,1)
%         %set stimulus template on the times
%         x = [x; t_start(t) + sig_stim_temp(:,1)];
%         y = [y; sig_stim_temp(:,2)]; %functionvalue
%         
%         x_cs = [x_cs; t_start(t) + sig_stim_temp(:,1)];
%         y_cs = [y_cs; sig_stim_temp(:,2)/cs_factor]; 
%         if ismember(t,l_id) %add atrial answer fitting to the s2 cv
%             cnts2 = cnts2 + 1;
%             if cnts2 < ERPid
%                 x = [x; t_start(t) + LATS2(el,cnts2) + atrial_sig_temp(:,1)];
%                 y = [y; atrial_sig_temp(:,2)];
%                 t_aa_s2(cnts2,el) = [t_start(t) + LATS2(el,cnts2)];
%                 
%                 x_cs = [x_cs; t_start(t) + LATS2(el,cnts2) + cs_window_delay + atrial_sig_temp(:,1)];
%                 y_cs = [y_cs; atrial_sig_temp(:,2)*3*cs_factor];
%             else
%                 t_aa_s2(cnts2,el) = nan;
%             end
%         else %add atrail answer fitting to s1 cv
%             x = [x; t_start(t) + LATS1(el) + atrial_sig_temp(:,1)];
%             y = [y; atrial_sig_temp(:,2)];
%             t_aa_s1(t-cnts2,el) = [t_start(t) + LATS1(el)];
%             
%             x_cs = [x_cs; t_start(t) + LATS1(el) + cs_window_delay + atrial_sig_temp(:,1)];
%             y_cs = [y_cs; atrial_sig_temp(:,2)*3*cs_factor]; 
%         end
%     end
%     
%     %add buffer to start & end of signal
%     x_buf(:,1) = linspace(0,buffer/1000,buffer/1000 * samplerate);
%     x = [x_buf ; x ; x(end) + 1/samplerate + x_buf];
%     y = [repmat(0,[buffer/1000*samplerate 1]) ; y ; repmat(0,[buffer/1000*samplerate 1]) ]; %buffervalues
%     signal{el} = [x y];
%     
%     x_cs = [x_buf ; x_cs ; x(end) + 1/samplerate + x_buf];
%     y_cs = [repmat(0,[buffer/1000*samplerate 1]) ; y_cs ; repmat(0,[buffer/1000*samplerate 1]) ];
%     signal_cs{el} = [x_cs y_cs];
% end
% 
% %set a last point so that all traces have same length
% for el = 1:size(LATS1,1)
%     signal{el}(end+1,1) = t_start(end) + x_buf(end) + s2s1dist;
%     signal{el}(end,2) = 0;
%     
%     signal_cs{el}(end+1,1) = t_start(end) + x_buf(end) + s2s1dist;
%     signal_cs{el}(end,2) = 0;
% end
% 
% %go through signal and if atrial answer slopes out before its stimulus cut
% %it off, else sum up signals -> overlap
% %probably creating empty signal zeros(size,1) and adding signals and AA is
% %better aproach
% [C,ia,ic] = unique(signal{1}(:,1));
% 
% 
% %signal(isnan(signal)) = 0;
% %add buffer to front & back of signal
% 
% 
% %for el = 1:size(LATS1,1)
% %figure
% %plot(signal{el}(:,1),signal{el}(:,2))
% %end
% %use interp if u need continous signal