classdef signal < handle
% -------------------------------------------------------
%
%    signal  - class for storing signal data
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
% Signal Class used as standard datastructure for signals throughout the program.
% By just inputting ydata & samplerate you can create a plot going from
% 1 to bla. If you want to have a reference time you can either input it
% as xdata or use the function setx1withsam
%
% setx1withsam sets the first x value as a value you define and
% calculates all the further values from the samplerate given.
%
% Inputs:
% ydata      - (timesteps x channels)
% samplerate - sampling of signal in (Hz) ,e.g 1000
%
% Outputs:
%           signal class
%
%
% Example Usage:
%    signal('name','mytestsignal','ydata',rand(500,3),'samplerate',500)
%   
% Revision history:

    
    
    properties
        ydata
        samplerate
        name
        id
        leads
        elecs
        yscale
        yax_name
        xax_name
        xdata
    end
    
    methods
        function obj = signal(varargin)
            p = inputParser;
            addParameter(p,'ydata','',@(x) ~isempty(x) && isnumeric(x) && numel(x)>1);
            addParameter(p,'samplerate','',@(x) ~isempty(x) && isnumeric(x) && numel(x)==1);
            addParameter(p,'name','',@(x) ~isempty(x) && ischar(x));
            addParameter(p,'leads','',@(x) ~isempty(x) && iscellstr(x));
            addParameter(p,'elecs','',@(x) ~isempty(x) && iscellstr(x));
            
            addParameter(p,'id','',@(x) ~isempty(x) && mod(x,1)==0);
            addParameter(p,'yax_name','',@(x) ~isempty(x) && ischar(x));
            addParameter(p,'xax_name','',@(x) ~isempty(x) && ischar(x));
            addParameter(p,'yscale','',@(x) ~isempty(x) && isnumeric(x)); %FOR NOW ONLY!!!! THIS HAS TO BE INCLUDED LATER!
            addParameter(p,'xdata','',@(x) ~isempty(x) && isnumeric(x)); %FOR NOW ONLY!!!! THIS HAS TO BE INCLUDED LATER!
            
            try
                p.parse(varargin{:});
            catch MyError
                rethrow(MyError);
            end
            
            obj.ydata = p.Results.ydata;
            obj.samplerate = p.Results.samplerate;
            obj.name = p.Results.name;
            if size(p.Results.leads,1)>size(p.Results.leads,2)
                obj.leads = p.Results.leads;
            else
                obj.leads = p.Results.leads';
            end
            if size(p.Results.elecs,1)>size(p.Results.elecs,2)
                obj.elecs = p.Results.elecs;
            else
                obj.elecs = p.Results.elecs';
            end
            obj.id = p.Results.id;
            obj.yax_name = p.Results.yax_name;
            obj.xax_name = p.Results.xax_name;
            obj.yscale = p.Results.yscale;
            obj.xdata = p.Results.xdata;
            
        end
        
        function setx1withsam(obj,xstart,samplerate)
            %setx1withsam sets the first x value as a value you define and
            %calculates all the further values from the samplerate given.
            if size(obj.xdata,1) > 1
                obj.xdata = [xstart:1/samplerate:size(obj.ydata,1)]';
            else
                fprintf('Signalclass: Imported xdata has less than 2 values. Import longer datastream.\n');
            end
        end
        
        function setsignal(obj,ydata, samplerate, name, id, leads, yax_name, xax_name, yscale, xdata)
            if nargin-1 == 0
                %fprintf('Signalclass: No input detected. Creating empty instance\n')
                obj.ydata = [];
                obj.samplerate = [];
                obj.name = [];
                obj.id = [];
                obj.leads = [];
                obj.yax_name = [];
                obj.xax_name = [];
                obj.yscale = [];
                obj.xdata = [];
                return
            elseif nargin-1 == 1
                fprintf('Signalclass: Only y-datapoints given and no samplingrate. Assuming samplerate of 1Hz!\n')
                if size(ydata,1)>size(ydata,2)
                    obj.ydata = ydata;
                else
                    obj.ydata = ydata';
                end
                obj.samplerate = 1;
                obj.name = [];
                obj.id = 0;
                obj.leads = [];
                obj.yax_name = [];
                obj.xax_name = [];
                obj.yscale = 1;
                obj.xdata = [1:size(ydata,2)]';
            elseif nargin-1 == 2
                if size(ydata,1)>size(ydata,2)
                    obj.ydata = ydata;
                else
                    obj.ydata = ydata';
                end
                obj.samplerate = samplerate;
                obj.name = [];
                obj.id = 0;
                obj.leads = [];
                obj.yax_name = [];
                obj.xax_name = [];
                obj.yscale = 1;
                obj.xdata = [(0:size(obj.ydata,1)-1)/samplerate]';
            elseif nargin-1 == 3
                if size(ydata,1)>size(ydata,2)
                    obj.ydata = ydata;
                else
                    obj.ydata = ydata';
                end
                obj.samplerate = samplerate;
                obj.name = name;
                obj.id = 0;
                obj.leads = [];
                obj.yax_name = [];
                obj.xax_name = [];
                obj.yscale = 1;
                obj.xdata = [(0:size(obj.ydata,1)-1)/samplerate]';
            elseif nargin-1 == 4
                if size(ydata,1)>size(ydata,2)
                    obj.ydata = ydata;
                else
                    obj.ydata = ydata';
                end
                obj.samplerate = samplerate;
                obj.name = name;
                obj.id = id;
                obj.leads = [];
                obj.yax_name = [];
                obj.xax_name = [];
                obj.yscale = 1;
                obj.xdata = [(0:size(obj.ydata,1)-1)/samplerate]';
            elseif nargin-1 == 5
                if size(ydata,1)>size(ydata,2)
                    obj.ydata = ydata;
                else
                    obj.ydata = ydata';
                end
                obj.samplerate = samplerate;
                obj.name = name;
                obj.id = id;
                obj.leads = leads;
                obj.yax_name = [];
                obj.xax_name = [];
                obj.yscale = 1;
                obj.xdata = [(0:size(obj.ydata,1)-1)/samplerate]';
            elseif nargin-1 == 6
                if size(ydata,1)>size(ydata,2)
                    obj.ydata = ydata;
                else
                    obj.ydata = ydata';
                end
                obj.samplerate = samplerate;
                obj.name = name;
                obj.id = id;
                obj.leads = leads;
                obj.yax_name = yax_name;
                obj.xax_name = [];
                obj.yscale = 1;
                obj.xdata = [(0:size(obj.ydata,1)-1)/samplerate]';
            elseif nargin-1 == 7
                if size(ydata,1)>size(ydata,2)
                    obj.ydata = ydata;
                else
                    obj.ydata = ydata';
                end
                obj.samplerate = samplerate;
                obj.name = name;
                obj.id = id;
                obj.leads = leads;
                obj.yax_name = yax_name;
                obj.xax_name = xax_name;
                obj.yscale = 1;
                obj.xdata = [(0:size(obj.ydata,1)-1)/samplerate]';
            elseif nargin-1 == 11
                if size(ydata,1)>size(ydata,2)
                    obj.ydata = ydata;
                else
                    obj.ydata = ydata';
                end
                obj.samplerate = samplerate;
                obj.name = name;
                obj.id = id;
                obj.leads = leads;
                obj.yax_name = yax_name;
                obj.xax_name = xax_name;
                obj.yscale = yscale;
                obj.xdata = [(0:size(obj.ydata,1)-1)/samplerate]';
            elseif nargin-1 == 12
                if size(ydata,1)>size(ydata,2)
                    obj.ydata = ydata;
                else
                    obj.ydata = ydata';
                end
                obj.samplerate = abs(xdata(1,:)-xdata(2,:));
                obj.name = name;
                obj.id = id;
                obj.leads = leads;
                obj.yax_name = yax_name;
                obj.xax_name = xax_name;
                obj.yscale = yscale;
                if size(xdata,1)>size(xdata,2)
                    obj.xdata = xdata;
                else
                    obj.xdata = xdata';
                end
            end
            
            
            
            
        end
    end
    
end

