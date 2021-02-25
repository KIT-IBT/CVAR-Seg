% -------------------------------------------------------
%
%    cretes unipolar templates based on dipole in infinite conductor method
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

function [pulse_final,t_shift] = createUnipolarAtrialSigTemplate_dipole(cv,x0,measuring_electrode,t_end,samplerate,amp)
%CREATEUNIPOLARATRIALSIGTEMPLATE_DIPOLE Summary of this function goes here
%   Detailed explanation goes here
%returned pulse is centered oround 0
% createUnipolarAtrialSigTemplate_dipole([500;0],[0;0],[20;0],0.200,samplerate,1);

cv = [cv; 0];            %mm/s
x0 = [x0; 0];
z_dist = 0.1;
measuring_electrode = [measuring_electrode;z_dist]; %adding z component sice elc is above tissue
cv_norm = norm(cv);
%t_new = 2*pdist2(x0',measuring_electrode')/cv_norm;

% t_end = 0.200;            %should be no longer than 200 ms
% samplerate = 2000;
% amp = 1;                  %mV;
edge_len    = [1; 1; 0];
elec_points = [measuring_electrode+edge_len/2 , measuring_electrode-edge_len/2 , measuring_electrode+[-1; +1; 1].*edge_len/2 , measuring_electrode+[+1; -1; 1].*edge_len/2];
%elec_points = measuring_electrode;

sigma       = sqrt(norm(cv));
dt          = 1/samplerate;
t           = 0:dt:t_end;%t_new;%   %s
y_fun       = @(cv,t,x0) [cv.*t + x0];
p           = repmat([1;1;0],[1 length(t)]);

for elp = 1:size(elec_points,2)
y           = y_fun(cv,t,x0);
dist_vec    = elec_points(:,elp)-y;
d           = pdist2(y',elec_points(:,elp)')';
u           = dist_vec/norm(dist_vec);

phi(elp,:)  = dot(p,u)./(4*pi*sigma*d.^2);%
end
phi = mean(phi);

%scale to desired max amplitude
scale_factor = amp/max(abs(phi));
phi          = phi * scale_factor;

%find timestep where dipole closest to measurementpoint & shift 2 zero
[~,t_in] = min(pdist2(y',measuring_electrode'));
t_shift = t'-t_in*dt;

% hold on
% plot(t*1000,phi,'-')
% plot(t*1000,mean(phi),'-')
%final curve
pulse_final  = [t', phi'];

%figure
%hold on
%plot(pulse_final(:,1),pulse_final(:,2))
end

