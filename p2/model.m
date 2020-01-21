close all
clear all
addpath('C:\Users\user\Desktop\code\control_lab\')
load('C:\Users\user\Desktop\code\control_lab\p2\M12_Grupo3_Presion_datos_depurados.mat')
path = 'C:\Users\user\Desktop\code\control_lab\p2\';
cd(path)

input = vop;
output = pv_smooth;

clearvars -except input offset_input output offset_output time
fig = 1;
% remove offset

input = input - offset_input;
output = output - offset_output;

% % make time uniform by interpolation
% mean sampling frequency was 2.2 samples per sec
time_uniform = min(time):1/2.2:max(time); 
output_interpolated = interp1(time,output,time_uniform);
input_interpolated = interp1(time,input,time_uniform);
error_interp = immse(output_interpolated,output(1:end-2));

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

% taking into account perturbation
[stableTime,stableVal,stableIndex] = findStablePoint(time,output,40,0.005);

% not taking into account perturbation
%[stableTime,stableVal,stableIndex] = findStablePoint(time,output,5,0.005);

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
output_zn = step(sys_zn,time)';



% Method: Miller Modification
tm_miller = tm_zn;
y = output(stableIndex)*(1-exp(-1));
tau_miller = (y - b_zn)/m_zn - tm_miller;
sys_miller = (kp/(tau_miller*s + 1))*exp(-tm_miller*s);
output_miller = step(sys_miller,time)';

% Method: Smith
[~,t28] = findClosest(output,output(stableIndex)*0.283); 
t28 = time(t28);
[~,t63] = findClosest(output,output(stableIndex)*0.632); 
t63 = time(t63);
tau_smith = 1.5*(t63-t28);
tm_smith = t63 - tau_smith;

sys_smith = (kp/(tau_smith*s + 1))*exp(-tm_smith*s);
output_smith = step(sys_smith,time)';

% Method: Alfaro
p1 = 0.25;
p2 = 0.75;
a = -0.91;
b = -1*a;
c = 1.262;
d = -1*(c - 1);
[tm_alfaro,tau_alfaro] = model2points(time,output,output(stableIndex),p1,p2,a,b,c,d);
sys_alfaro = (kp/(tau_smith*s + 1))*exp(-tm_smith*s);
output_alfaro = step(sys_alfaro,time)';

% Method: Broida
p1 = 0.28;
p2 = 0.40;
a = -5.5;
b = -1*a;
c = 2.8;
d = -1*(c - 1);
[tm_broida,tau_broida] = model2points(time,output,output(stableIndex),p1,p2,a,b,c,d);
sys_broida = (kp/(tau_broida*s + 1))*exp(-tm_broida*s);
output_broida = step(sys_broida,time)';

% Method Chen - Yang
p1 = 0.33;
p2 = 0.67;
a = -1.4;
b = -1*a;
c = 1.540;
d = -1*(c - 1);
[tm_cy,tau_cy] = model2points(time,output,output(stableIndex),p1,p2,a,b,c,d);
sys_cy = (kp/(tau_cy*s + 1))*exp(-tm_cy*s);
output_cy = step(sys_cy,time)';

% Method Ho et al

p1 = 0.35;
p2 = 0.85;
a = -0.670;
b = -1*a;
c = 1.300;
d = -0.290;
[tm_ho,tau_ho] = model2points(time,output,output(stableIndex),p1,p2,a,b,c,d);
sys_ho = (kp/(tau_ho*s + 1))*exp(-tm_ho*s);
output_ho = step(sys_ho,time)';

% Method Viteckova et al

p1 = 0.33;
p2 = 0.70;
a = -1.245;
b = -1*a;
c = 1.498;
d = -1*(c - 1);
[tm_viteckova,tau_viteckova] = model2points(time,output,output(stableIndex),p1,p2,a,b,c,d);
sys_viteckova = (kp/(tau_viteckova*s + 1))*exp(-tm_viteckova*s);
output_viteckova = step(sys_viteckova,time)';

% plot responses
fig = fig + 1;
figure(fig)
plot(time,output_zn,'r-','LineWidth',2);
hold on
plot(time,output_miller,'g-','LineWidth',2);
hold on
plot(time,output_smith,'b-','LineWidth',2);
hold on
plot(time,output_alfaro,'y-','LineWidth',2);
hold on
plot(time,output_broida,'c-','LineWidth',2);
hold on
plot(time,output_cy,'m-','LineWidth',2);
hold on
plot(time,output_ho,'r:','LineWidth',2);
hold on
plot(time,output_viteckova,'b:','LineWidth',2);
hold on

% we need to rescale the original because step was not unitary
% we could use time_raw and output_raw here if we want...
plot(time_raw,output_raw/max(input),'k-','LineWidth',2);
title('Comparación modelos 1er orden')
legend('Ziegler-Nichols','Miller','Smith','Alfaro','Broida','Chen-Yang','Ho et al','Viteckova et al','Real','Location','SE')
xlabel('Tiempo (s)')
saveas(gcf,'comparison_1st_order','png');

% MSE for 1st order models
model_names_1 = {'Alfaro','Broida','Chen-Yang','Ho et al','Miller','Smith','Viteckova et al','Ziegler-Nichols'};
output_models_1= [output_alfaro;output_broida;output_cy;output_ho;output_miller;output_smith;output_viteckova;output_zn];
numberOfModels_1 = 8;
MSE_1 = zeros(8,1);
for i = 1:numberOfModels_1
    MSE_1(i) = immse(output_models_1(i,:),output);
end

[error_1,idx_1] = min(MSE_1);
best_1 = model_names_1(idx_1);


