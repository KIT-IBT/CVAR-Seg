% -------------------------------------------------------
%
%    plotMultiChannelEGM  - Code to plot multichannel EGMs
%
%    Ver. 0.1.0 (alpha)
%
%    Created:       Tobias Oesterlein        (06.12.2015)
%    Last modified: Tobias Oesterlein        (06.12.2015)
%
%    Institute of Biomedical Engineering
%    Universitaet Karlsruhe (TH)
%
%    http://www.ibt.uni-karlsruhe.de
%
%    Based on filtering techniques from Gustavo Lenis.
%
%    Copyright 2000-2015 - All rights reserved.
%
% ------------------------------------------------------
%
% [fh] = plotMultiChannelEGM(sigplot)
% Simple way to plot multichannel EMG data
%
% Input: structure sigplot with fields
%        [title]:           String being figure title
%
%        [fs]:           	sampling rate in Hz. If provided, sample axis
%                               will be replaced by axis in seconds
%
%        [samplespan]:      Samples of the signal to be plotted
%
%        [position]:        Position vector for figure (4 x 1)
%
%        [LineWidth]:       LineWidth for figure (scalar)
%
%        sig:               Cell array containg each signal in fields
%         [leads]            names of N leads (N x 1 cell array)
%
%         [prefix]           name prefix for this catheter (1 x 1 char array)
%
%         signals            signal data of N leads (T samples x N)
%
%         [ampFactor]        amplitfication for these leads after normalization (scalar)
%
%         [color]            color of N leads (character or 3x1 array)
%
%         [WFM]              Wavefront matrix with indices of W wavefronts (matrix NxW)
%
%
% Output:
%        fh:              handle to figure window
%
%
% Example Usage:
%   sigplot.title='MiFi signals during Ablation';
% 
%   sigplot.sig{1}.leads={'I' 'II' 'V6'};
%   sigplot.sig{1}.signals=ECG;
%   sigplot.sig{1}.color='k';
% 
%   sigplot.sig{2}.leads={'ABL12' 'ABL34'};
%   sigplot.sig{2}.signals=signals(:,[1 3]);
%   sigplot.sig{2}.color=[0 0 0.8];
%  
%   plotMultiChannelEGM(sigplot);
% 
%
% Revision:
%
%
%


function [fh,normfact] = plotMultiChannelEGM(sigplot)

%% check input variables

if ~isfield(sigplot,'position')
    positiondata=[100 100 400 400];
else
    positiondata=sigplot.position;
end
    
if ~isfield(sigplot,'LineWidth')
    LineWidthdata=1;
else
    LineWidthdata=sigplot.LineWidth;
end

if ~isfield(sigplot,'title')
    titledata='Channel plot';
else
    titledata=sigplot.title;
end

if ~isfield(sigplot,'fs') || isempty(sigplot.fs)
    samplingRateGiven=false;
else
    samplingRateGiven=true;
    samplingRate=sigplot.fs;
end

Nt=size(sigplot.sig{1}.signals,1);
if ~isfield(sigplot,'samplespan') || isempty(sigplot.samplespan)
    samplespandata=1:Nt;
else
    samplespandata=sigplot.samplespan;
end

Ns=length(sigplot.sig);
listDefaultColors=prism(Ns);
listDefaultColors(ismember(listDefaultColors,[1 1 0],'rows'),:)=0;    % replace yellow by black
for nSig=1:Ns
    % make sure that all signals have correct dimensions
    Nl=size(sigplot.sig{nSig}.signals,2);
    assert(size(sigplot.sig{nSig}.signals,1)==Nt,...
        sprintf('ERROR: Number of samples in signal %d is not %d',nSig,Nt));    
    
    % set values to default if required: color
    if ~isfield(sigplot.sig{nSig},'color') || isempty(sigplot.sig{nSig}.color)
        sigplot.sig{nSig}.color=listDefaultColors(nSig,:);
    end
    
    % set values to default if required: amplitude scaling factor
    if ~isfield(sigplot.sig{nSig},'ampFactor') || isempty(sigplot.sig{nSig}.ampFactor)
        sigplot.sig{nSig}.ampFactor=1;
    end
    
    % set values to default: Name of leads for labeling
    if ~isfield(sigplot.sig{nSig},'leads') || isempty(sigplot.sig{nSig}.leads)
        sigplot.sig{nSig}.leads=cell(Nl,1);
        for n1=1:Nl
            sigplot.sig{nSig}.leads{n1}=sprintf('C %d',n1);
        end
    else
        % entries in field leads are present, so check them
        if length(sigplot.sig{nSig}.leads)==Nl
            % make sure that dimension is correct for vertical adding
            sigplot.sig{nSig}.leads=reshape(sigplot.sig{nSig}.leads,[],1);
            % add the prefix if it is present
            if isfield(sigplot.sig{nSig},'prefix') && ~isempty(sigplot.sig{nSig}.prefix)
                for nLead=1:Nl
                    sigplot.sig{nSig}.leads(nLead)={[sigplot.sig{nSig}.prefix sigplot.sig{nSig}.leads{nLead}]};
                end
                clear nLead
            end
        elseif length(sigplot.sig{nSig}.leads)==1
            % copy the first entry and number if
            sigprefix=sigplot.sig{nSig}.leads{1};
            sigplot.sig{nSig}.leads=cell(Nl,1);
            for n1=1:Nl
                sigplot.sig{nSig}.leads{n1}=sprintf([sigprefix ' %d'],n1);
            end
        else
             error('ERROR: Number of traces in signal %d does not match number of leads',nSig);
        end            
    end
    
    % make sure that number of channels is correct
    if ~isfield(sigplot.sig{nSig},'WFM')
        sigplot.sig{nSig}.WFM=[];
    end
    if ~isempty(sigplot.sig{nSig}.WFM)
        if size(sigplot.sig{nSig}.WFM,1)~=Nl
            warning('Number of channels in WFM of signal %d does not match. Skipping entry ...',nSig)
            sigplot.sig{nSig}.WFM=[];
        end
    end
end

% time vector
if samplingRateGiven
    timevec=(0:Nt-1)/samplingRate;
else
    timevec=1:Nt;
end
  


%% generate plot

% initialize plot
fh=figure('Position',positiondata);
legendlist=[];
offsetCount=0;

% add one signal group by each other
for nSig=1:Ns
    % get data
    sig=sigplot.sig{nSig};
    
    % compute offset
    datain=sig.signals(samplespandata,:);
    numSigs=size(datain,2);
    legendlist=vertcat(legendlist,strrep(sig.leads,'_',' '));
    offset=repmat(offsetCount+(0:numSigs-1),size(datain,1),1);
    normfact=max(abs(datain));
    normfact(isnan(normfact))=1;
    normfact(normfact==0)=1;
    datanorm=datain./repmat(normfact,size(datain,1),1)*sig.ampFactor;
    
    %plot(timevec(samplespandata),datain-offset,'Color',sig.color,'LineWidth',LineWidthdata)
    plot(timevec(samplespandata),datanorm-offset,'Color',sig.color,'LineWidth',LineWidthdata)
    offsetCount=offsetCount+numSigs;
    hold on
    
    % add wavefron if provided
    if ~isempty(sig.WFM)
        for nWF=1:size(sig.WFM,2)
            % wavefronts
            sigWithLAT=find(~isnan(sig.WFM(:,nWF)));
            LATinds=sig.WFM(sigWithLAT,nWF);
            LATtimes=timevec(LATinds);
            % individual LATs
            plot(LATtimes,-offset(1,sigWithLAT),'xr');
            % connected LATs
            tmpPattern=zeros(size(offset(1,:)));
            tmpPattern(sigWithLAT)=1;
            tmpCC=bwconncomp(tmpPattern);
            for nCC=1:tmpCC.NumObjects
                if length(tmpCC.PixelIdxList{nCC})<2; continue; end
                plot(timevec(sig.WFM(tmpCC.PixelIdxList{nCC},nWF)),-offset(1,tmpCC.PixelIdxList{nCC}),'-r');
            end
        end
    end
end

% finalize plot
grid on
title(titledata)
set(gca,'YTick',[-(offsetCount-1):0],'YTickLabel',flipud(legendlist))

if samplingRateGiven
    xlabel('Time [seconds]');
else
    xlabel('Time [samples]');
end

axis tight

ylim([-offsetCount-2 3])





