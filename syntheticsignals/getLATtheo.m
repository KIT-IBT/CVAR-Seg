% -------------------------------------------------------
%
%    getLATtheo - gets theretical LATs based on elliptic propagation
%
%    Ver. 1.0.0
%
%    Created:           Mark Nothstein, Johannes Tischer (25.02.2020)
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

function [LAT,dist,cath_mes,vec_ang] = getLATtheo(def_cat_rad,cv_long,cv_trans,elips_rot)

% location of bipolar meassurement/stimulation in x- and y-direction
cath_mes= nan(2,10);
for c = 0:9
    cath_mes(:,c+1) = [def_cat_rad*cos((18+c*36)*pi/180)  def_cat_rad*sin((18+c*36)*pi/180)];
end
% distance between bipolar stimulation and bipolar meassurement
%dist = nan(1,30);
for c = 1:10
    dist(c,1) = sqrt((cath_mes(1,c) - cath_mes(1,2))^2 + (cath_mes(2,c) - cath_mes(2,2))^2);
    %dist(1,c+10) = sqrt((cath_mes(1,c) - cath_mes(1,6))^2 + (cath_mes(2,c) - cath_mes(2,6))^2);
    %dist(1,c+20) = sqrt((cath_mes(1,c) - cath_mes(1,9))^2 + (cath_mes(2,c) - cath_mes(2,9))^2);
end
% calculate the direction of the vetor in rad
vec_ang = zeros(10,1);
for c = 1:10
    if (cath_mes(1,c) - cath_mes(1,2)) > 0
        vec_ang(c,1) = atan((cath_mes(2,c)-cath_mes(2,2))/(cath_mes(1,c)-cath_mes(1,2)));
    elseif (cath_mes(1,c) - cath_mes(1,2)) < 0 && (cath_mes(2,c)-cath_mes(2,2)) >= 0
        vec_ang(c,1) = atan((cath_mes(2,c)-cath_mes(2,2))/(cath_mes(1,c)-cath_mes(1,2))) + pi ;
    elseif (cath_mes(1,c) - cath_mes(1,2)) < 0 && (cath_mes(2,c)-cath_mes(2,2)) < 0
        vec_ang(c,1) = atan((cath_mes(2,c)-cath_mes(2,2))/(cath_mes(1,c)-cath_mes(1,2))) - pi ;
    elseif (cath_mes(1,c) - cath_mes(1,2)) == 0 && (cath_mes(2,c)-cath_mes(2,2)) > 0
        vec_ang(c,1) = pi/2 ;
    elseif (cath_mes(1,c) - cath_mes(1,2)) == 0 && (cath_mes(2,c)-cath_mes(2,2)) < 0
        vec_ang(c,1) = - pi/2 ;
    end
end
cath_mes = cath_mes';
%
% for c = 1:10
%     if (cath_mes(1,c) - cath_mes(1,6)) > 0
%         vec_ang(c,2) = atan((cath_mes(2,c)-cath_mes(2,6))/(cath_mes(1,c)-cath_mes(1,6)));
%     elseif (cath_mes(1,c) - cath_mes(1,6)) < 0 && (cath_mes(2,c)-cath_mes(2,6)) >= 0
%         vec_ang(c,2) = atan((cath_mes(2,c)-cath_mes(2,6))/(cath_mes(1,c)-cath_mes(1,6))) + pi ;
%     elseif (cath_mes(1,c) - cath_mes(1,6)) < 0 && (cath_mes(2,c)-cath_mes(2,6)) < 0
%         vec_ang(c,2) = atan((cath_mes(2,c)-cath_mes(2,6))/(cath_mes(1,c)-cath_mes(1,6))) - pi ;
%     elseif (cath_mes(1,c) - cath_mes(1,6)) == 0 && (cath_mes(2,c)-cath_mes(2,6)) > 0
%         vec_ang(c,2) = pi/2 ;
%     elseif (cath_mes(1,c) - cath_mes(1,6)) == 0 && (cath_mes(2,c)-cath_mes(2,6)) < 0
%         vec_ang(c,2) = - pi/2 ;
%     end
% end
%
% for c = 1:10
%     if (cath_mes(1,c) - cath_mes(1,9)) > 0
%         vec_ang(c,3) = atan((cath_mes(2,c)-cath_mes(2,9))/(cath_mes(1,c)-cath_mes(1,9)));
%     elseif (cath_mes(1,c) - cath_mes(1,9)) < 0 && (cath_mes(2,c)-cath_mes(2,9)) >= 0
%         vec_ang(c,3) = atan((cath_mes(2,c)-cath_mes(2,9))/(cath_mes(1,c)-cath_mes(1,9))) + pi ;
%     elseif (cath_mes(1,c) - cath_mes(1,9)) < 0 && (cath_mes(2,c)-cath_mes(2,9)) < 0
%         vec_ang(c,3) = atan((cath_mes(2,c)-cath_mes(2,9))/(cath_mes(1,c)-cath_mes(1,9))) - pi ;
%     elseif (cath_mes(1,c) - cath_mes(1,9)) == 0 && (cath_mes(2,c)-cath_mes(2,9)) > 0
%         vec_ang(c,3) = pi/2 ;
%     elseif (cath_mes(1,c) - cath_mes(1,9)) == 0 && (cath_mes(2,c)-cath_mes(2,9)) < 0
%         vec_ang(c,3) = - pi/2 ;
%     end
% end

%vec_ang = [vec_ang(:,1).', vec_ang(:,2).', vec_ang(:,3).'].' ;
%vec_ang = rmmissing(vec_ang);
%dist = dist(dist ~= 0 );

LAT = nan(size(vec_ang,1),1);
a = nan(size(vec_ang,1),1);
b = nan(size(vec_ang,1),1);
% calculation of the local activation time
for c = 1:size(vec_ang,1)
    a(c,1) = abs(dist(c,1) .* cos(vec_ang(c,1)));
    b(c,1) = abs(dist(c,1) .* sin(vec_ang(c,1)));
    LAT(c,1) = sqrt((a(c,1)/cv_long).^2 + (b(c,1)/cv_trans).^2) ;
end

end