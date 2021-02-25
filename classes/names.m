% -------------------------------------------------------
%
%    names  - class for viabl catheter names
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
%
% Revision history:

classdef names
    %CATHNAMES lookuptable for allowed catheter names of the program
    
    properties
        cathnames
    end
    
    methods
        function obj = names
            %   Add additional catheternames here
            obj.cathnames = ...
                {'CS';...
                'ABL';...
                'SPI';...
                'MAP';...
                'TEMPERATURE'
                'Lasso'...
                };
        end
        
        function logic_check = ispartof(obj,name_in)
            %checks input agains class strings
            if ischar(name_in) && size(name_in,1)==1
                logic_check = sum(strcmp(obj.cathnames,name_in));
            else
                fprintf('Wrong namestring. Check naming conventions in names class.\n');
            end
            
        end
    end
end

