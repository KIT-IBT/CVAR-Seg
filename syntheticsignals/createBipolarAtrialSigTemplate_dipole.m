% -------------------------------------------------------
%
%    createBipolarAtrialSigTemplate_dipole - creates synthethic bipolar atrial
%    activity from 2 unipolar signals 
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

function pulse_final = createBipolarAtrialSigTemplate_dipole(cv,x0,measuring_electrode_1,measuring_electrode_2,t_end,samplerate,amp)
%
% createBipolarAtrialSigTemplate_dipole([500;0],[0;0],[3;0],[20;0],0.200,samplerate,1)

%show path
figure
hold on
plot(x0(1),x0(2),'*k')
plot(measuring_electrode_1(1),measuring_electrode_1(2),'*r')
plot(measuring_electrode_2(1),measuring_electrode_2(2),'*r')
t = [0: 1/(samplerate) :t_end ];
res = cv*t(end);
quiver(x0(1),x0(2),res(1),res(2));
xlim([min([x0(1) measuring_electrode_1(1) measuring_electrode_2(1)]) max([x0(1) measuring_electrode_1(1) measuring_electrode_2(1)])]);
%ylim([min([x0(2) measuring_electrode_1(2) measuring_electrode_2(2)]) max([x0(2) measuring_electrode_1(2) measuring_electrode_2(2)])]);
%plot(t,y_line(t),'-r')
hold off

[el_1, t_shift_1] = createUnipolarAtrialSigTemplate_dipole(cv,x0,measuring_electrode_1,t_end,samplerate,amp);
[el_2, t_shift_2] = createUnipolarAtrialSigTemplate_dipole(cv,x0,measuring_electrode_2,t_end,samplerate,amp);

bip_sig = el_1(:,2)-el_2(:,2);
%scale to desired max amplitude
scale_factor = amp/max(abs(bip_sig));

pulse_final = [el_1(:,1),bip_sig*scale_factor];

figure
plot(pulse_final(:,1),pulse_final(:,2))
end