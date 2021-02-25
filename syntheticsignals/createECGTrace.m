% -------------------------------------------------------
%
%    createECGtrace - creates ECG trace - ToDo
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

function ECGTrace = createECGTrace(t_start,sigTrace)
%uses precreated ECG Signal and sets beginning of p-wave to t_start
    ECGTrace(:,1)=sigTrace;
    ECGTrace(:,2)=zeros(size(sigTrace,1),1);
end
