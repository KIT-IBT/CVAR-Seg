classdef timesegment < handle
    %SEGMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        min_fs
        max_fs
        
        peaks
        numofpeaks
        peaks_bip_fs
        numofpeaks_bip
        peaks_uni_fs
        numofpeaks_uni
        peaktype        %unipolar or bipolar -> how was it found
        
        content
        
        segment_id
        segtype         %What kind of peak -> S1 , S2 or nothing
        segtype_num
        
        weight
        
        channelid_postsplit     %needed to reconstruct where everysegment originated
        
        history         %cell of strings explaining what was done e.g HPfilt,NLEO,
    end
    
    properties (Access = private)
        
    end
    
    methods
        function obj = timesegment(minsample,maxsample,peak,varargin)
            %Constructs
            p = inputParser;
            addRequired(p,'minsample',@(x) ~isempty(x) && isnumeric(x));
            addRequired(p,'maxsample',@(x) ~isempty(x) && isnumeric(x));
            addRequired(p,'peak',@(x) ~isempty(x) && isnumeric(x));
            addParameter(p,'content','',@(x) ~isempty(x) && (isnumeric(x) || isa(x,'timesegment')) );
            addParameter(p,'weight','',@(x) ~isempty(x) && isnumeric(x));
            addParameter(p,'chanid_postsplit','',@(x) ~isempty(x) && isnumeric(x) && size(x,1) == 1 && size(x,2) == 1);
            addParameter(p,'peaktype','',@(x) ischar(x) && (strcmp(x,'unipolar') || strcmp(x,'bipolar')));
            try
                p.parse( minsample,maxsample,peak,varargin{:} );
            catch MyError
                rethrow(MyError);
            end
            
            obj.min_fs = p.Results.minsample;
            obj.max_fs = p.Results.maxsample;
            obj.peaks = p.Results.peak;
            obj.content = p.Results.content;
            obj.weight = p.Results.weight;
            obj.channelid_postsplit = p.Results.chanid_postsplit;
            %obj.peaktype = p.Results.peaktype;
            %if strcmp(obj.peaktype,'unipolar')
            %    obj.peaks_uni_fs = p.Results.peak;
            %elseif strcmp(obj.peaktype,'bipolar')
            %    obj.peaks_bip_fs = p.Results.peak;
            %end
            
        
            
            
        end
        
        function addpeak(obj,peak_fs)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.updatenumofpeaks;
            obj.peaks(obj.numofpeaks+1) = peak_fs;
        end
        
        function addpeak_bipolar(obj,peak_fs)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.peaks_bip_fs = peak_fs;
        end
        
        function updatenumofpeaks(obj)
            obj.numofpeaks = size(obj.peaks,1);
        end
        
        function addcontent(obj,content_in)
            if isa(content_in,'timesegment')
                if isempty(obj.content)
                    obj.content{1} = content_in;
                else
                    obj.content{end+1} = content_in;
                end
            else
                fprintf('wrong inputformat.\n');
            end
        end
        
        function addsegtype(obj,string)
            if ischar(string)
            obj.segtype = string;
            else
                fprintf('Input is not string \n')
            end
        end
        
    end
end

