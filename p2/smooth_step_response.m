% Data was preprocessed in Excel prior to importing to matlab
% So as to generate the time vector accurately 
% given that the sampling frequency is not constant
% Basically each second has its own sampling frequency
clear all
% Loading of the data
path = 'C:\Users\user\Desktop\code\control_lab\p2\';
filename = 'to_matlab.csv';
data = csvread([path filename]);

% Partitioning
time_all = data(:,1)';  % time vector
pv_all = data(:,2)'; % process value, in this case pressure
% set value is column 3, not needed since there is no control being done
vop_all = data(:,4)'; % valve opening percentage

% find when the step stimuli starts
start = find(vop_all>=max(vop_all), 1, 'first'); 
time = time_all(start:end);
time = time - time_all(start); % remove offset to start from zero
pv = pv_all(start:end);
vop = vop_all(start:end);

% Smoothing
span = 0.1 * length(time); % span of 10% of the total number of data points
pv_smooth = smooth(time,pv,span,'moving')'; % smoothing


% Plotting
figure(1)
plot(time,pv,'g-','LineWidth',3)
hold on 
plot(time,pv_smooth,'b-')
title('Respuesta Experimental vs Suavizada')
legend('Respuesta de presión cruda','Respuesta de presión suavizada','Location','SE')
xlabel('Tiempo (s)')
ylabel('Presión ante estimulo de tipo escalón en t=0 (psig)')

save([path 'M12_Grupo3_Presion_datos_depurados.mat'])
csvwrite([path 'M12_Grupo3_Presion_datos_depurados.csv'],[time' vop' pv' pv_smooth'])