% Data was preprocessed in Excel prior to importing to matlab
% So as to generate the time vector accurately 
% given that the sampling frequency is not constant
% Basically each second has its own sampling frequency
clear all
% Loading of the data
path = 'Y:\y\sync\current\courses\control\lab\p2\';
filename = 'to_matlab.csv';
data = csvread([path filename]);

% Partitioning
time_all = data(:,1)';  % time vector
pv_all = data(:,2)'; % process value, in this case pressure
% set value is column 3, not needed since there is no control being done
vop_all = data(:,4)'; % valve opening percentage


start = find(vop_all>=max(vop_all), 1, 'first');
time = time_all(start:end);
time = time - time_all(start); % remove offset to start from zero
pv = pv_all(start:end);
vop = vop_all(start:end);


% Normalization for easier visualization
pv_norm = pv/norm(pv);
vop_norm = vop/norm(vop);

% Smoothing
span = 0.1 * length(time); % span of 5% of the total number of data points
pv_smooth = smooth(time,pv,span,'moving'); % smoothing
pv_smooth_norm = pv_smooth/norm(pv_smooth);

% Plotting
figure(1)
plot(time,vop_norm,'r-','LineWidth',2)
hold on
plot(time,pv_norm,'g-','LineWidth',3)
hold on 
plot(time,pv_smooth_norm,'b-')

legend('Estimulo de apertura de valvula','Respuesta de presi�n cruda','Respuesta de presi�n suavizada','Location','SE')
xlabel('Tiempo (s)')
ylabel('Valores normalizados para comparaci�n')

figure(2)
plot(time,pv,'g-','LineWidth',3)
hold on 
plot(time,pv_smooth,'b-')

legend('Respuesta de presi�n cruda','Respuesta de presi�n suavizada','Location','SE')
xlabel('Tiempo (s)')
ylabel('Presi�n')

save('M12_Grupo3_Presion_datos_depurados.mat','pv_smooth')