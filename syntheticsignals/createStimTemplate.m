% -------------------------------------------------------
%
%    createStimTemplate - creates synthethic stimulation template
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

function pulse_final = createStimTemplate(d1,D,S,s2,tau_edge,tau_plateau,B,samplerate)

dt = 1/samplerate; %ms samplerate
buffer_factor_biphase_shift = 5;
buffer_factor_endpoint = 5;


t1 = d1*D - buffer_factor_biphase_shift * tau_edge;
t2 = D - buffer_factor_endpoint * tau_edge;
tf1 = linspace(0,t1,1/dt);
tf2 = linspace(t1+dt,t2,1/dt);
tf3 = linspace(t2+dt,D,1/dt);
tfend_t = D+dt;
t = [tf1 tf2 tf3 tfend_t] - t1;

S2=-s2*S;
f1 = B + S .* (1-exp( -tf1./tau_edge )) .* exp( -tf1./tau_plateau );
f2 = B + S2 .* (1-exp( -(tf2-t1)./tau_edge )) .* exp( -(tf2-t1)./tau_plateau );
Pt2_exact = B + S2 .* (1-exp( -(t2-t1)./tau_edge )) .* exp( -(t2-t1)./tau_plateau );
f3 = B + Pt2_exact * exp( -(tf3-t2)./tau_edge );
fend = B;
f_fin = [f1 f2 f3 fend]';

x_uniform = linspace(tf1(1,1),tfend_t,abs(tf1(1,1)-tfend_t)/1000 * samplerate)' - t1;
f_uniform(:,1) = interp1(t,f_fin,x_uniform);

%final curve
pulse_final = [t' f_fin];
end