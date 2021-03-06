% Data was preprocessed in Excel prior to importing to matlab
% So as to generate the time vector accurately 
% given that the sampling frequency is not constant
% Basically each second has its own sampling frequency
close all
clear all
addpath('y:\code\control_lab\')
% Loading of the data
path = 'Y:\code\control_lab\p2\data\';
cd(path)
filename = 'to_matlab.csv';
data = csvread(filename);

% Partitioning
time_all = data(:,1)';  % time vector
pv_all = data(:,2)'; % process value, in this case pressure
% set value is column 3, not needed since there is no control being done
vop_all = data(:,4)'; % valve opening percentage


% map to current

a = 4; b = 19.2;

vop_all = mapScale(vop_all,0,100,a,b);

[stableTime,stableVal,stableIndex] = findStablePoint(time_all,pv_all,40,0.005);
tm = time_all(find(pv_all>pv_all(stableIndex), 1, 'first')); 


% find when the step stimuli starts
start = find(vop_all>=max(vop_all), 1, 'first'); 
tm = tm - time_all(start);
time = time_all(start:end);
time = time - time_all(start); % remove offset to start from zero
pv = pv_all(start:end);
vop = vop_all(start:end);



% Smoothing
span = 0.1 * length(time); % span of 10% of the total number of data points
pv_smooth = smooth(time,pv,span,'moving')'; % smoothing


% Plotting
fig = 1;
figure(fig)
plot(time,pv,'g-','LineWidth',3)
hold on 
plot(time,pv_smooth,'k-','LineWidth',2)
title('Respuesta Experimental vs Suavizada')
legend('Respuesta de presi�n cruda','Respuesta de presi�n suavizada','Location','SE')
xlabel('Tiempo (s)')
ylabel('Presi�n ante estimulo de tipo escal�n en t=0 (mA)')
grid on
saveas(gcf,'smooth_vs_experimental','png');

fig = fig + 1;
figure(fig)
plot(time_all,pv_all,'b-','LineWidth',2)
hold on 
plot(time_all,vop_all,'r-','LineWidth',2)
title('Porcentaje de Apertura de la Valvula (entrada) vs Presi�n del Tanque (salida)')
legend('Salida (mA)','Entrada (mA)','Location','SE')
xlabel('Tiempo (s)')
grid on
saveas(gcf,'input_vs_output_raw','png');

fig = fig + 1;
figure(fig) % Normalized input/output plot
plot(time_all,pv_all/norm(pv_all),'b-','LineWidth',2)
hold on 
plot(time_all,vop_all/norm(vop_all),'r-','LineWidth',2)
title('Porcentaje de Apertura de la Valvula (entrada) vs Presi�n del Tanque (salida) normalizadas')
legend('Salida Normalizada','Entrada Normalizada','Location','SE')
xlabel('Tiempo (s)')
grid on
saveas(gcf,'input_vs_output_raw_norm','png');

fig = fig + 1;
figure(fig) % Normalized input/output plot
plot(time,pv_smooth/norm(pv_smooth),'b-','LineWidth',2)
hold on 
plot(time,vop/norm(vop),'r-','LineWidth',2)
title('Porcentaje de Apertura de la Valvula (entrada) vs Presi�n del Tanque (salida) normalizadas y suavizada')
legend('Salida Normalizada Suavizada','Entrada Normalizada','Location','SE')
xlabel('Tiempo (s)')
grid on
saveas(gcf,'input_vs_output_norm_smooth','png');

offset_input = min(vop_all);
offset_output = min(pv_smooth); % important that the offset uses the signal used after, not the original

pv_smooth_all = [min(pv_smooth)*ones(1,length(pv_all)-length(pv_smooth)),pv_smooth];

save([path 'M12_Grupo3_Presion_datos_depurados.mat'])
csvwrite([path 'M12_Grupo3_Presion_datos_depurados.csv'],[time' vop' pv' pv_smooth'])