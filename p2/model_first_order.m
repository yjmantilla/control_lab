close all
clear all
addpath('C:\Users\user\Desktop\code\control_lab\')
load('C:\Users\user\Desktop\code\control_lab\p2\M12_Grupo3_Presion_datos_depurados.mat')
path = 'C:\Users\user\Desktop\code\control_lab\p2\';
cd(path)

input = vop;
output = pv_smooth;

clearvars -except input offset_input output offset_output time

% remove offset

input = input - offset_input;
output = output - offset_output;

% make time uniform by interpolation
% mean sampling frequency was 2.2 samples per sec
time_uniform = min(time):1/2.2:max(time); 
output_interpolated = interp1(time,output,time_uniform);
input_interpolated = interp1(time,input,time_uniform);

fig =1;
figure(fig)
plot(time,output/norm(output),'g-','LineWidth',4)
hold on
plot(time_uniform,output_interpolated/norm(output_interpolated),'k-','LineWidth',2)
%hold on
%plot(time_uniform,input_interpolated/norm(input_interpolated),'k-','LineWidth',2)
title('Generación de tiempo uniforme por interpolación')
legend('Salida no uniforme','Salida Interpolada Uniforme','Location','SE')
xlabel('Tiempo (s)')
saveas(gcf,'interpolation','png');

time_raw = time;
output_raw = output;
input_raw = input;
time = time_uniform;
output = output_interpolated;
input = input_interpolated;
% Input Output Plot

fig = fig + 1;
figure(fig)
plot(time,input,'b-','LineWidth',2)
hold on 
plot(time,output,'r-','LineWidth',2)
title('Entrada vs Salida sin estado estable')
legend('Entrada (% de apertura valvula)','Salida (presión psig)','Location','SE')
xlabel('Tiempo (s)')
saveas(gcf,'input_vs_output_no_offset','png');

fig = fig + 1;
figure(fig)
plot(time,input/norm(input),'b-','LineWidth',2)
hold on 
plot(time,output/norm(output),'r-','LineWidth',2)
title('Entrada vs Salida sin estado estable, normalizadas')
legend('Entrada (% de apertura valvula normalizada)','Salida (presión psig normalizada)','Location','SE')
xlabel('Tiempo (s)')
saveas(gcf,'input_vs_output_no_offset_norm','png');

% find when the signal is stable and determine the value
[stableTime,stableVal,stableIndex] = findStablePoint(time,output,40);

% find gain (notice that offset was removed so no need to do a delta)

kp = output(stableIndex)/input(stableIndex);

% Method: Ziegler Nichols
% since time is not equispaced, better to use diff
slope = diff(output)./diff(time);
[maxSlope,maxSlope_index] = max(slope);

fig = fig + 1;
figure(fig)
plot(time(1:end-1),slope)
title('Derivada numerica (diff(y)./diff(t)) de la respuesta respecto al tiempo (no uniforme)')

m_zn = maxSlope;
b_zn = output(maxSlope_index) - m_zn*time(maxSlope_index);
%tangent = maxSlope*(time-time(maxSlope_index))+output(maxSlope_index);
tangent = m_zn * time + b_zn;

% tm : time between stimuli start and cut of the x axis by the tangent
y = 0;
tm_zn = (y-b_zn)/m_zn;
% tau : time between tm and where the tangent reaches the final value of
% the response
y = output(stableIndex);
% tau_zn =(time(maxSlope_index) + (max(output)-output(maxSlope_index))/maxSlope)-tm_zn; 
tau_zn = (y-b_zn)/m_zn - tm_zn;

fig = fig + 1;
figure(fig)
plot(time,tangent,'k-','LineWidth',2);
hold on
plot(time,output,'r-','LineWidth',2);
ylim([min(output) max(output)])
title('Recta de mayor tangente vs Respuesta')
legend('Recta Tangente','Respuesta','Location','SE')
xlabel('Tiempo (s)')
saveas(gcf,'zn_tangent','png');

% Create the system
s = tf('s');
sys_zn = (kp/(tau_zn*s + 1))*exp(-tm_zn*s);
output_zn = step(sys_zn,time);



% Method: Miller Modification
tm_miller = tm_zn;
y = output(stableIndex)*(1-exp(-1));
tau_miller = (y - b_zn)/m_zn - tm_miller;
sys_miller = (kp/(tau_miller*s + 1))*exp(-tm_miller*s);
output_miller = step(sys_miller,time);


% plot responses
fig = fig + 1;
figure(fig)
plot(time,output_zn,'r-','LineWidth',2);
hold on
plot(time,output_miller,'g-','LineWidth',2);
hold on
% we need to rescale the original because step was not unitary
plot(time,output/max(input),'k-','LineWidth',2);
title('Comparación modelos 1er orden')
legend('Ziegler Nichols','Miller','Real','Location','SE')
xlabel('Tiempo (s)')
saveas(gcf,'comparison_1st_order','png');