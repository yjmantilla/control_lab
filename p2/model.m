close all
clear all
addpath('C:\Users\user\Desktop\code\control_lab\')
load('C:\Users\user\Desktop\code\control_lab\p2\M12_Grupo3_Presion_datos_depurados.mat')
path = 'C:\Users\user\Desktop\code\control_lab\p2\';
cd(path)

plant_in = vop;
plant_out = pv_smooth;

clearvars -except plant_in offset_input plant_out offset_output time
aux_fig = 1;
% remove offset

plant_in = plant_in - offset_input;
plant_out = plant_out - offset_output;

% % make time uniform by interpolation
% mean sampling frequency was 2.2 samples per sec
interpolated_time = min(time):1/2.2:max(time); 
interpolated_out = interp1(time,plant_out,interpolated_time);
interpolated_in = interp1(time,plant_in,interpolated_time);
interpolated_error = immse(interpolated_out,plant_out(1:end-2));

figure(aux_fig)
plot(time,plant_out/norm(plant_out),'g-','LineWidth',4)
hold on
plot(interpolated_time,interpolated_out/norm(interpolated_out),'k-','LineWidth',2)
%hold on
%plot(time_uniform,in_interpolated/norm(in_interpolated),'k-','LineWidth',2)
title('Generaci�n de tiempo uniforme por interpolaci�n')
legend('Salida no uniforme','Salida Interpolada Uniforme','Location','SE')
xlabel('Tiempo (s)')
saveas(gcf,'interpolation','png');

raw_time = time;
raw_out = plant_out;
raw_in = plant_in;
time = interpolated_time;
plant_out = interpolated_out;
plant_in = interpolated_in;


% in out Plot

aux_fig = aux_fig + 1;
figure(aux_fig)
plot(time,plant_in,'b-','LineWidth',2)
hold on 
plot(time,plant_out,'r-','LineWidth',2)
title('Entrada vs Salida sin estado estable')
legend('Entrada (% de apertura valvula)','Salida (presi�n psig)','Location','SE')
xlabel('Tiempo (s)')
saveas(gcf,'in_vs_out_no_offset','png');

aux_fig = aux_fig + 1;
figure(aux_fig)
plot(time,plant_in/norm(plant_in),'b-','LineWidth',2)
hold on 
plot(time,plant_out/norm(plant_out),'r-','LineWidth',2)
title('Entrada vs Salida sin estado estable, normalizadas')
legend('Entrada (% de apertura valvula normalizada)','Salida (presi�n psig normalizada)','Location','SE')
xlabel('Tiempo (s)')
saveas(gcf,'in_vs_out_no_offset_norm','png');

% find when the signal is stable and determine the value

% taking into account perturbation
[plant_stableTime,plant_stableVal,plant_stableIndex] = findStablePoint(time,plant_out,40,0.005);

% not taking into account perturbation
%[plant_stableTime,plant_stableVal,plant_stableIndex] = findStablePoint(time,out,5,0.005);

% find gain (notice that offset was removed so no need to do a delta)

plant_kp = plant_out(plant_stableIndex)/plant_in(plant_stableIndex);

% 1st Order Methods

% Method: Ziegler Nichols
% since time is not equispaced, better to use diff
zn_slope = diff(plant_out)./diff(time);
[zn_maxSlope,zn_maxSlope_index] = max(zn_slope);

aux_fig = aux_fig + 1;
figure(aux_fig)
plot(time(1:end-1),zn_slope)
title('Derivada numerica (diff(y)./diff(t)) de la respuesta respecto al tiempo (no uniforme)')

zn_m = zn_maxSlope;
zn_b = plant_out(zn_maxSlope_index) - zn_m*time(zn_maxSlope_index);
%tangent = maxSlope*(time-time(maxSlope_index))+out(maxSlope_index);
zn_tangent = zn_m * time + zn_b;

% tm : time between stimuli start and cut of the x axis by the tangent
aux_y = 0;
zn_tm = (aux_y-zn_b)/zn_m;
% tau : time between tm and where the tangent reaches the final value of
% the response
aux_y = plant_out(plant_stableIndex);
% tau_zn =(time(maxSlope_index) + (max(out)-out(maxSlope_index))/maxSlope)-tm_zn; 
zn_tau = (aux_y-zn_b)/zn_m - zn_tm;

aux_fig = aux_fig + 1;
figure(aux_fig)
plot(time,zn_tangent,'k-','LineWidth',2);
hold on
plot(time,plant_out,'r-','LineWidth',2);
ylim([min(plant_out) max(plant_out)])
title('Recta de mayor tangente vs Respuesta')
legend('Recta Tangente','Respuesta','Location','SE')
xlabel('Tiempo (s)')
saveas(gcf,'zn_tangent','png');

% Create the system
s = tf('s');
zn_sys = (plant_kp/(zn_tau*s + 1))*exp(-zn_tm*s);
zn_out = step(zn_sys,time)';



% Method: Miller Modification
miller_tm = zn_tm;
aux_y = plant_out(plant_stableIndex)*(1-exp(-1));
miller_tau = (aux_y - zn_b)/zn_m - miller_tm;
miller_sys = (plant_kp/(miller_tau*s + 1))*exp(-miller_tm*s);
miller_out = step(miller_sys,time)';

% Method: Smith
[~,t28] = findClosest(plant_out,plant_out(plant_stableIndex)*0.283); 
t28 = time(t28);
[~,t63] = findClosest(plant_out,plant_out(plant_stableIndex)*0.632); 
t63 = time(t63);
smith_tau = 1.5*(t63-t28);
smith_tm = t63 - smith_tau;

smith_sys = (plant_kp/(smith_tau*s + 1))*exp(-smith_tm*s);
smith_out = step(smith_sys,time)';

% Method: Alfaro
m2p_p1 = 0.25;
m2p_p2 = 0.75;
m2p_a = -0.91;
m2p_b = -1*m2p_a;
m2p_c = 1.262;
m2p_d = -1*(m2p_c - 1);
[alfaro_tm,alfaro_tau] = model2points(time,plant_out,plant_out(plant_stableIndex),m2p_p1,m2p_p2,m2p_a,m2p_b,m2p_c,m2p_d);
alfaro_sys = (plant_kp/(alfaro_tau*s + 1))*exp(-alfaro_tm*s);
alfaro_out = step(alfaro_sys,time)';

% Method: Broida
m2p_p1 = 0.28;
m2p_p2 = 0.40;
m2p_a = -5.5;
m2p_b = -1*m2p_a;
m2p_c = 2.8;
m2p_d = -1*(m2p_c - 1);
[broida_tm,broida_tau] = model2points(time,plant_out,plant_out(plant_stableIndex),m2p_p1,m2p_p2,m2p_a,m2p_b,m2p_c,m2p_d);
broida_sys = (plant_kp/(broida_tau*s + 1))*exp(-broida_tm*s);
broida_out = step(broida_sys,time)';

% Method Chen - Yang
m2p_p1 = 0.33;
m2p_p2 = 0.67;
m2p_a = -1.4;
m2p_b = -1*m2p_a;
m2p_c = 1.540;
m2p_d = -1*(m2p_c - 1);
[cy_tm,cy_tau] = model2points(time,plant_out,plant_out(plant_stableIndex),m2p_p1,m2p_p2,m2p_a,m2p_b,m2p_c,m2p_d);
cy_sys = (plant_kp/(cy_tau*s + 1))*exp(-cy_tm*s);
cy_out = step(cy_sys,time)';

% Method Ho et al

m2p_p1 = 0.35;
m2p_p2 = 0.85;
m2p_a = -0.670;
m2p_b = -1*m2p_a;
m2p_c = 1.300;
m2p_d = -0.290;
[ho_tm,ho_tau] = model2points(time,plant_out,plant_out(plant_stableIndex),m2p_p1,m2p_p2,m2p_a,m2p_b,m2p_c,m2p_d);
ho_sys = (plant_kp/(ho_tau*s + 1))*exp(-ho_tm*s);
ho_out = step(ho_sys,time)';

% Method Viteckova et al

m2p_p1 = 0.33;
m2p_p2 = 0.70;
m2p_a = -1.245;
m2p_b = -1*m2p_a;
m2p_c = 1.498;
m2p_d = -1*(m2p_c - 1);
[viteckova_tm,viteckova_tau] = model2points(time,plant_out,plant_out(plant_stableIndex),m2p_p1,m2p_p2,m2p_a,m2p_b,m2p_c,m2p_d);
viteckova_sys = (plant_kp/(viteckova_tau*s + 1))*exp(-viteckova_tm*s);
viteckova_out = step(viteckova_sys,time)';

% plot responses
aux_fig = aux_fig + 1;
figure(aux_fig)
plot(time,zn_out,'r-','LineWidth',2);
hold on
plot(time,miller_out,'g-','LineWidth',2);
hold on
plot(time,smith_out,'b-','LineWidth',2);
hold on
plot(time,alfaro_out,'y-','LineWidth',2);
hold on
plot(time,broida_out,'c-','LineWidth',2);
hold on
plot(time,cy_out,'m-','LineWidth',2);
hold on
plot(time,ho_out,'r:','LineWidth',2);
hold on
plot(time,viteckova_out,'b:','LineWidth',2);
hold on

% we need to rescale the original because step was not unitary
% we could use time_raw and out_raw here if we want...
plot(time,plant_out/max(plant_in),'k-','LineWidth',2);
title('Comparaci�n modelos 1er orden')
legend('Ziegler-Nichols','Miller','Smith','Alfaro 1/4 3/4','Broida','Chen-Yang','Ho et al','Viteckova et al','Real','Location','SE')
xlabel('Tiempo (s)')
saveas(gcf,'comparison_1st_order','png');

% MSE for 1st order models
models_1_names = {'Alfaro 1/4 3/4','Broida','Chen-Yang','Ho et al','Miller','Smith','Viteckova et al','Ziegler-Nichols'};
models_1_outs= [alfaro_out;broida_out;cy_out;ho_out;miller_out;smith_out;viteckova_out;zn_out];
models_1_number = 8;
models_1_MSE = zeros(8,1);
%remember to compare against scaled out for unit step
for aux_i = 1:models_1_number
    models_1_MSE(aux_i) = immse(models_1_outs(aux_i,:),plant_out/max(plant_in));
end

[models_1_error,models_1_idx] = min(models_1_MSE);
models_1_best = models_1_names(models_1_idx);


% 2nd Order Methods

% Method: Alfaro General 123c
[~,t25] = findClosest(plant_out,plant_out(plant_stableIndex)*0.25); 
t25 = time(t25);
[~,t50] = findClosest(plant_out,plant_out(plant_stableIndex)*0.5); 
t50 = time(t50);
[~,t75] = findClosest(plant_out,plant_out(plant_stableIndex)*0.75); 
t75 = time(t75);
g123c_a = 23;
g123c_tauPP = (t75-t25)/(0.9866+0.7036*g123c_a);
%tau1_123c = 