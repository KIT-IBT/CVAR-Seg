% -------------------------------------------------------
%
%    ismultiple - check if values is a multiple of the other within error
%    range
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

function [ismultiple,multiplier_exact] = ismultiple(y1,y2,err)

ymax = max([y1,y2]);
ymin = min([y1,y2]);
largestval = ceil(ymax/ymin);

remainder = 0;
ismultiple = 0;
multiplier_exact = 0;

n = 1;
for i=1:largestval
    %remainder = mod(ymax/(i*ymin),1);
    remainder = abs(ymax-(i*ymin));
    if remainder < err
        multiplier_exact = ymax/(ymin);
        ismultiple = 1; 
    end
end

end
