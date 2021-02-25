% -------------------------------------------------------
%
%    add_baselinewander - ads synthetic baselinewander
%     
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

function [sigTrace_final] = add_baselinewander(sigTrace,samplerate,baseline_mean_amp,for_each_chan)
%sigTrace_final must be samples x channels

    power = baseline_mean_amp; %amp
    x_uniform = (0:1/samplerate:(length(sigTrace)-1)/samplerate);
    vals = [0.001 0.005 0.01 0.02 0.03 0.05 0.1 0.25 0.5 0.75 1 1.25 1.5 1.75 2 2.5 3 3.5 4 4.5 5];
    bline = zeros(size(sigTrace,1),2,2);
    for v = 1:numel(vals)
        if for_each_chan ~=0
            for chan = 1:size(sigTrace,2)
                phi = rand(1);
                p_rand(v,chan) = power * 2 * rand(1); %amp
                bline(:,v,chan) = p_rand(v,chan) * cos(2*pi*vals(v)*x_uniform/1000 + 2*pi*phi);
                sigTrace_final(:,chan) = sigTrace(:,chan) + bline(:,v,chan);
            end
        else
            p_rand(v) = power * 2 * rand(1); %amp
            phi = rand(1);
            bline(:,v,1) = p_rand(v) * cos(2*pi*vals(v)*x_uniform/1000 + 2*pi*phi);
            sigTrace_final = sigTrace + bline(:,v,1);
        end
    end
    
end

