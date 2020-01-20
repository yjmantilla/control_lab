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
% Input Output Plot

figure(1)
plot(time,input,'b-','LineWidth',2)
hold on 
plot(time,output,'r-','LineWidth',2)
title('Entrada vs Salida sin estado estable')
legend('Entrada (% de apertura valvula)','Salida (presión psig)','Location','SE')
xlabel('Tiempo (s)')
saveas(gcf,'input_vs_output_no_offset','png');

figure(2)
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

% Ziegler Nichols
% since time is not equispaced, better to use diff
slope = diff(output)./diff(time);
[maxSlope,maxSlope_index] = max(slope);
figure(3)
plot(time(1:end-1),slope)
title('Derivada numerica (diff(y)./diff(t)) de la respuesta respecto al tiempo (no uniforme)')
tangent = maxSlope*(time-time(maxSlope_index))+output(maxSlope_index);
tm_zn = time(maxSlope_index) - output(maxSlope_index)/maxSlope;
tau_zn = (time(maxSlope_index) + (max(output)-output(maxSlope_index))/maxSlope)-tm_zn; 
figure(4)
plot(time,tangent,'k-','LineWidth',2);
hold on
plot(time,output,'r-','LineWidth',2);
ylim([min(output) max(output)])
title('Recta de mayor tangente vs Respuesta')
legend('Recta Tangente','Respuesta','Location','SE')
xlabel('Tiempo (s)')
saveas(gcf,'zn_tangent','png');