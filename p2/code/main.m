close all
clear all
addpath('Y:\code\control_lab\')
load('Y:\code\control_lab\p2\data\M12_Grupo3_Presion_datos_depurados.mat')
path = 'Y:\code\control_lab\p2\';
cd(path)

plant.in = vop;
plant.out = pv_smooth;
plant.in_all = vop_all;
plant.out_all = pv_smooth_all;
time.t = time;
time.all = time_all;
offset.input = offset_input;
offset.output = offset_output;
clearvars -except plant offset time
aux_fig = 1;

% remove offset

plant.in = plant.in - offset.input;
plant.out = plant.out - offset.output;
plant.in_all = plant.in_all - offset.input;
plant.out_all = plant.out_all - offset.output;

% % make time uniform by interpolation
% mean sampling frequency was 2.2 samples per sec
time.fs = 2.2;
time.Ts = 1/time.fs;
time.ending = 400;
interpolated.time = min(time.t):time.Ts:max(time.t);
interpolated.out = interp1(time.t,plant.out,interpolated.time);
interpolated.in = interp1(time.t,plant.in,interpolated.time);
interpolated.error = immse(interpolated.out,plant.out(1:end-2));


interpolated.time_all = min(time.all):time.Ts:max(time.all);
interpolated.out_all = interp1(time.all,plant.out_all,interpolated.time_all);
interpolated.in_all = interp1(time.all,plant.in_all,interpolated.time_all);

interpolated.time_ext = min(time.all):time.Ts:time.ending;
figure(aux_fig)
plot(time.t,plant.out/norm(plant.out),'g-','LineWidth',4)
hold on
plot(interpolated.time,interpolated.out/norm(interpolated.out),'k-','LineWidth',2)
%hold on
%plot(time.uniform,in_interpolated/norm(in_interpolated),'k-','LineWidth',2)
title('Generación de tiempo uniforme por interpolación')
legend('Salida no uniforme','Salida Interpolada Uniforme','Location','SE')
xlabel('Tiempo (s)')
grid on
saveas(gcf,'interpolation','png');

raw.time = time.t;
raw.out = plant.out;
raw.in = plant.in;
time.t = interpolated.time;
plant.out = interpolated.out;
plant.in = interpolated.in;

time.all = interpolated.time_all;
plant.out_all = interpolated.out_all;
plant.in_all = interpolated.in_all;

% Extension of plant response assuming it keeps it stable value
time.ext = interpolated.time_ext;
plant.out_ext = [plant.out ones(1,length(time.ext)-length(plant.out))*max(plant.out)];

% in out Plot

aux_fig = aux_fig + 1;
figure(aux_fig)
plot(time.t,plant.in,'b-','LineWidth',2)
hold on 
plot(time.t,plant.out,'r-','LineWidth',2)
title('Entrada vs Salida sin estado estable')
legend('Entrada (mA apertura valvula)','Salida (presión mA)','Location','SE')
xlabel('Tiempo (s)')
grid on
saveas(gcf,'in_vs_out_no_offset','png');

aux_fig = aux_fig + 1;
figure(aux_fig)
plot(time.t,plant.in/norm(plant.in),'b-','LineWidth',2)
hold on 
plot(time.t,plant.out/norm(plant.out),'r-','LineWidth',2)
title('Entrada vs Salida sin estado estable, normalizadas')
legend('Entrada (mA apertura valvula normalizada)','Salida (presión mA normalizada)','Location','SE')
xlabel('Tiempo (s)')
grid on
saveas(gcf,'in_vs_out_no_offset_norm','png');

% find when the signal is stable and determine the value

% taking into account perturbation
[plant.stableTime,plant.stableVal,plant.stableIndex] = findStablePoint(time.t,plant.out,40,0.005);

% not taking into account perturbation
%[plant.stableTime,plant.stableVal,plant.stableIndex] = findStablePoint(time.t,plant.out,5,0.005);

% find gain (notice that offset was removed so no need to do a delta)

plant.kp = plant.out(plant.stableIndex)/plant.in(plant.stableIndex);

% 1st Order Methods

% Method: Ziegler Nichols
% since time is not equispaced, better to use diff
zn.name = 'Ziegler - Nichols';
zn.slope = diff(plant.out)./diff(time.t);
zn.slope = gradient(plant.out,time.Ts);
[zn.maxSlope,zn.maxSlope_index] = max(zn.slope);

aux_fig = aux_fig + 1;
figure(aux_fig)
plot(time.t(1:length(zn.slope)),zn.slope)
grid on
title('Derivada numerica (diff(y)./diff(t)) de la respuesta respecto al tiempo (no uniforme)')

zn.m = zn.maxSlope;
zn.b = plant.out(zn.maxSlope_index) - zn.m*time.t(zn.maxSlope_index);
%tangent = maxSlope*(time.t-time.t(maxSlope_index))+out(maxSlope_index);
zn.tangent = zn.m * time.t + zn.b;

% tm : time between stimuli start and cut of the x axis by the tangent
zn.y = 0;
zn.tm = (zn.y-zn.b)/zn.m;
% tau : time between tm and where the tangent reaches the final value of
% the response
zn.y = plant.out(plant.stableIndex);
% tau_zn =(time.t(maxSlope_index) + (max(out)-out(maxSlope_index))/maxSlope)-tm_zn; 
zn.tau = (zn.y-zn.b)/zn.m - zn.tm;

aux_fig = aux_fig + 1;
figure(aux_fig)
plot(time.t,zn.tangent,'k-','LineWidth',2);
hold on
plot(time.t,plant.out,'r-','LineWidth',2);
ylim([min(plant.out) max(plant.out)])
title('Recta de mayor tangente vs Respuesta')
legend('Recta Tangente','Respuesta','Location','SE')
xlabel('Tiempo (s)')
grid on
saveas(gcf,'zn.tangent','png');

% Create the system
s = tf('s');
zn.sys = (plant.kp/(zn.tau*s + 1))*exp(-zn.tm*s);
zn.out = step(zn.sys,time.ext)';



% Method: Miller Modification
miller.name = 'Miller';
miller.tm = zn.tm;
miller.y = plant.out(plant.stableIndex)*(1-exp(-1));
%miller.tau = (miller.y - zn.b)/zn.m - miller.tm;
[~,miller.tau] = findClosest(plant.out,miller.y);
miller.tau = time.t(miller.tau)-miller.tm;
miller.sys = (plant.kp/(miller.tau*s + 1))*exp(-miller.tm*s);
miller.out = step(miller.sys,time.ext)';

% Method: Smith
smith.name = 'Smith';
[~,time.t28] = findClosest(plant.out,plant.out(plant.stableIndex)*0.283); 
time.t28 = time.t(time.t28);
[~,time.t63] = findClosest(plant.out,plant.out(plant.stableIndex)*0.632); 
time.t63 = time.t(time.t63);
smith.tau = 1.5*(time.t63-time.t28);
smith.tm = time.t63 - smith.tau;

smith.sys = (plant.kp/(smith.tau*s + 1))*exp(-smith.tm*s);
smith.out = step(smith.sys,time.ext)';

% Method: Alfaro
alfaro.name = 'Alfaro 1/4 3/4';
m2p.p1 = 0.25;
m2p.p2 = 0.75;
m2p.a = -0.91;
m2p.b = -1*m2p.a;
m2p.c = 1.262;
m2p.d = -1*(m2p.c - 1);
[alfaro.tm,alfaro.tau] = model2points(time.t,plant.out,plant.out(plant.stableIndex),m2p.p1,m2p.p2,m2p.a,m2p.b,m2p.c,m2p.d);
alfaro.sys = (plant.kp/(alfaro.tau*s + 1))*exp(-alfaro.tm*s);
alfaro.out = step(alfaro.sys,time.ext)';

% Method: Broida
broida.name = 'Broida';
m2p.p1 = 0.28;
m2p.p2 = 0.40;
m2p.a = -5.5;
m2p.b = -1*m2p.a;
m2p.c = 2.8;
m2p.d = -1*(m2p.c - 1);
[broida.tm,broida.tau] = model2points(time.t,plant.out,plant.out(plant.stableIndex),m2p.p1,m2p.p2,m2p.a,m2p.b,m2p.c,m2p.d);
broida.sys = (plant.kp/(broida.tau*s + 1))*exp(-broida.tm*s);
broida.out = step(broida.sys,time.ext)';

% Method Chen - Yang
cy.name = 'Chen - Yang';
m2p.p1 = 0.33;
m2p.p2 = 0.67;
m2p.a = -1.4;
m2p.b = -1*m2p.a;
m2p.c = 1.540;
m2p.d = -1*(m2p.c - 1);
[cy.tm,cy.tau] = model2points(time.t,plant.out,plant.out(plant.stableIndex),m2p.p1,m2p.p2,m2p.a,m2p.b,m2p.c,m2p.d);
cy.sys = (plant.kp/(cy.tau*s + 1))*exp(-cy.tm*s);
cy.out = step(cy.sys,time.ext)';

% Method Ho et al
ho.name = 'Ho et al';
m2p.p1 = 0.35;
m2p.p2 = 0.85;
m2p.a = -0.670;
m2p.b = -1*m2p.a;
m2p.c = 1.300;
m2p.d = -0.290;
[ho.tm,ho.tau] = model2points(time.t,plant.out,plant.out(plant.stableIndex),m2p.p1,m2p.p2,m2p.a,m2p.b,m2p.c,m2p.d);
ho.sys = (plant.kp/(ho.tau*s + 1))*exp(-ho.tm*s);
ho.out = step(ho.sys,time.ext)';

% Method Viteckova et al
viteckova.name = 'Viteckova et al';
m2p.p1 = 0.33;
m2p.p2 = 0.70;
m2p.a = -1.245;
m2p.b = -1*m2p.a;
m2p.c = 1.498;
m2p.d = -1*(m2p.c - 1);
[viteckova.tm,viteckova.tau] = model2points(time.t,plant.out,plant.out(plant.stableIndex),m2p.p1,m2p.p2,m2p.a,m2p.b,m2p.c,m2p.d);
viteckova.sys = (plant.kp/(viteckova.tau*s + 1))*exp(-viteckova.tm*s);
viteckova.out = step(viteckova.sys,time.ext)';

% Matlab Identification 1st Order
id1.name = 'Matlab 1er Orden';
id1.data = iddata(plant.out_all',plant.in_all',time.Ts);
%id1.data = iddata(plant.out',plant.in',time.Ts);
% sys = tfest(data,np,nz,iodelay); if NaN estimate
opt = procestOptions('InitialCondition','zero');
id1.sys = procest(id1.data,'P1D',opt);
id1.sys.Kp = plant.kp;
%id2.sys.num = plant.kp;
%id1.tau = get_taus(id1.sys.den);
%for aux_i = 1:length(id1.tau)
%id1.sys.num = id1.sys.num / id1.tau(aux_i);
%end

id1.out = step(id1.sys,time.ext)';

%id1.kp = id1.sys.num;
%for aux_i = 1:length(id1.tau)
%id1.kp = id1.kp * id1.tau(aux_i);
%end
%id1.tm = id1.sys.ioDelay;


% plot responses
aux_fig = aux_fig + 1;
figure(aux_fig)
plot(time.ext,zn.out,'y-','LineWidth',2);
hold on
plot(time.ext,miller.out,'g-','LineWidth',2);
hold on
plot(time.ext,smith.out,'b-','LineWidth',2);
hold on
plot(time.ext,alfaro.out,'r-','LineWidth',2);
hold on
%plot(time.ext,broida.out,'r:','LineWidth',2);
%hold on
%plot(time.ext,cy.out,'g:','LineWidth',2);
%hold on
plot(time.ext,ho.out,'c-','LineWidth',2);
hold on
%plot(time.ext,viteckova.out,'b:','LineWidth',2);
%hold on
plot(time.ext,id1.out,'m-','LineWidth',2);
hold on

% we need to rescale the original because step was not unitary
% we could use time.raw and out_raw here if we want...
plot(time.ext,plant.out_ext/max(plant.in),'k-','LineWidth',2);
title('Comparación modelos 1er orden')
%legend(zn.name,miller.name,smith.name,alfaro.name,broida.name,cy.name,ho.name,viteckova.name,id1.name,'Real','Location','SE')
legend(zn.name,miller.name,smith.name,alfaro.name,ho.name,id1.name,'Real','Location','SE')
xlabel('Tiempo (s)')
grid on
saveas(gcf,'comparison_1st_order','png');

% MSE for 1st order models
models_1.names = {alfaro.name,broida.name,cy.name,ho.name,miller.name,smith.name,viteckova.name,zn.name,id1.name};
models_1.outs= [alfaro.out;broida.out;cy.out;ho.out;miller.out;smith.out;viteckova.out;zn.out;id1.out];
models_1.number = 9;
models_1.MSE = zeros(models_1.number,1);
%remember to compare against scaled out for unit step
for aux_i = 1:models_1.number
    models_1.MSE(aux_i) = immse(models_1.outs(aux_i,:),plant.out_ext/max(plant.in));
end

[models_1.error,models_1.idx] = min(models_1.MSE);
models_1.best = models_1.names(models_1.idx);


% 2nd Order Methods

% Method: Alfaro General 123c
g123c.name = 'Alfaro General 123c';
[~,time.t25] = findClosest(plant.out,plant.out(plant.stableIndex)*0.25); 
time.t25 = time.t(time.t25);
[~,time.t50] = findClosest(plant.out,plant.out(plant.stableIndex)*0.5); 
time.t50 = time.t(time.t50);
[~,time.t75] = findClosest(plant.out,plant.out(plant.stableIndex)*0.75); 
time.t75 = time.t(time.t75);
g123c.a = (-0.6240*time.t25 + 0.9866*time.t50 - 0.3626*time.t75)/(0.3533*time.t25 - 0.7036*time.t50 + 0.3503*time.t75);
g123c.a = abs(g123c.a);
g123c.tauPP = (time.t75-time.t25)/(0.9866+0.7036*g123c.a);
%g123c.tauPP = abs(g123c.tauPP);
if g123c.tauPP < 0
    g123c.tauPP = 0;
end
g123c.tau1 = g123c.tauPP;
g123c.tau2 = g123c.a*g123c.tauPP;
if g123c.tau1 < 0
    g123c.tau1 = 0;
end

if g123c.tau2 < 0
    g123c.tau2 = 0;
end
g123c.tm = time.t75 - (1.3421+1.3455*g123c.a)*g123c.tauPP;
%g123c.tm = abs(g123c.tm);
if g123c.tm < 0
    g123c.tm = 0;
end
g123c.sys = (plant.kp * exp(-g123c.tm*s))/((g123c.tau1*s+1)*(g123c.tau2*s+1));
g123c.out = step(g123c.sys,time.ext)';

% Method : Simetrico
sym.name = 'Simetrico';
sym.numberOfTests = 200;
[sym.tau1,sym.tau2,sym.tm,sym.x] = symmetricModel(sym.numberOfTests,plant.kp,plant.in,plant.out,plant.stableIndex,time.t,1e-3);
sym.sys = (plant.kp * exp(-sym.tm*s))/((sym.tau1*s+1)*(sym.tau2*s+1));
sym.out = step(sym.sys,time.ext)';

% Method : Harriott
har.name = 'Harriott';
har.tm = 0;
[~,time.t73] = findClosest(plant.out,plant.out(plant.stableIndex)*0.73); 
har.tausum = (time.t73-har.tm)/1.3;
har.t = har.tm + 0.5*har.tausum;

[~,har.t_index] = findClosest(time.t,har.t); 
har.ratio = plant.out(har.t_index)/plant.out(plant.stableIndex);
% RATIO IS GREATER THAN 0.4 SO THE HARRIOTT METHOD IS NOT APPLICABLE

% Method: Stark
[~,time.t15] = findClosest(plant.out,plant.out(plant.stableIndex)*0.15); 
time.t15 = time.t(time.t15);
[~,time.t45] = findClosest(plant.out,plant.out(plant.stableIndex)*0.45); 
time.t45 = time.t(time.t45);
[~,time.t75] = findClosest(plant.out,plant.out(plant.stableIndex)*0.75); 
time.t75 = time.t(time.t75);
stark.name = 'Stark';
stark.x = (time.t45-time.t15)/(time.t75-time.t15);
stark.damp = (0.0805-5.547*(0.475-stark.x)^2)/(stark.x-0.356);
if stark.damp <= 1
    stark.f2 = 0.708*(2.811)^stark.damp;
else
    stark.f2 = 2.6*stark.damp -0.6;
end

stark.wn = stark.f2/(time.t75-time.t15);
stark.f3 = 0.922 * 1.66^(stark.damp);
stark.tm = time.t45 - stark.f3/stark.wn;
%stark.tm = abs(stark.tm);
if stark.tm < 0
    stark.tm = 0;
end
stark.tau1 = (stark.damp + sqrt(stark.damp^2 - 1))/(stark.wn);
stark.tau2 = (stark.damp - sqrt(stark.damp^2 - 1))/(stark.wn);

stark.sys = (plant.kp * exp(-stark.tm*s))/((stark.tau1*s+1)*(stark.tau2*s+1));
stark.out = step(stark.sys,time.ext)';

% Method: Jahanmiri - Fallahi
jf.name = 'Jahanmiri - Fallahi';
[~,time.t2] = findClosest(plant.out,plant.out(plant.stableIndex)*0.02); 
time.t2 = time.t(time.t2);
[~,time.t5] = findClosest(plant.out,plant.out(plant.stableIndex)*0.05); 
time.t5 = time.t(time.t5);
[~,time.t90] = findClosest(plant.out,plant.out(plant.stableIndex)*0.9); 
time.t90 = time.t(time.t90);
[~,time.t70] = findClosest(plant.out,plant.out(plant.stableIndex)*0.7); 
time.t70 = time.t(time.t70);
jf.tms = [time.t2,time.t5];
jf.errors = zeros(length(jf.tms),1);

for aux_i = 1:length(jf.tms)
    [jf.out,~] = jf_model(time.t,plant.kp,jf.tms(aux_i),time.t70,time.t90);
    jf.errors(aux_i) = immse(jf.out,plant.out/max(plant.in));
end

[~,jf.index] = min(jf.errors);
jf.tm = jf.tms(jf.index);
[jf.out,jf.sys,jf.tau] = jf_model(time.ext,plant.kp,jf.tm,time.t70,time.t90);

% Method: Ho et al - Polo Doble
ho2.name = 'Ho et al - Polo Doble';
m2p.p1 = 0.35;
m2p.p2 = 0.85;
m2p.a = -0.463;
m2p.b = -1*m2p.a;
m2p.c = 1.574;
m2p.d = -1*(m2p.c - 1);
[ho2.tm,ho2.tau] = model2points(time.t,plant.out,plant.out(plant.stableIndex),m2p.p1,m2p.p2,m2p.a,m2p.b,m2p.c,m2p.d);
ho2.sys = (plant.kp/((ho2.tau*s + 1)^2))*exp(-ho2.tm*s);
ho2.out = step(ho2.sys,time.ext)';

% Method: Viteckova et al - Polo Doble

viteckova2.name = 'Viteckova et al - Polo Doble';
m2p.p1 = 0.33;
m2p.p2 = 0.70;
m2p.a = -0.749;
m2p.b = -1*m2p.a;
m2p.c = 1.937;
m2p.d = -1*(m2p.c - 1);
[viteckova2.tm,viteckova2.tau] = model2points(time.t,plant.out,plant.out(plant.stableIndex),m2p.p1,m2p.p2,m2p.a,m2p.b,m2p.c,m2p.d);
viteckova2.sys = (plant.kp/((viteckova2.tau*s + 1)^2))*exp(-viteckova2.tm*s);
viteckova2.out = step(viteckova2.sys,time.ext)';

% Matlab Identification 2nd Order
id2.name = 'Matlab 2do Orden';
opt = procestOptions('InitialCondition','zero');
id2.data = iddata(plant.out_all',plant.in_all',time.Ts);
%id2.data = iddata(plant.out',plant.in',time.Ts);
%sys = tfest(data,np,nz,iodelay); if NaN estimate
id2.sys = procest(id2.data,'P2D',opt);
id2.sys.Kp = plant.kp;
%id2.sys.num = plant.kp;
%id2.taus = get_taus(id2.sys.den);
%for aux_i = 1:length(id2.taus)
%id2.sys.num = id2.sys.num / id2.taus(aux_i);
%end

%id2.kp = id2.sys.num;
%for aux_i = 1:length(id2.taus)
%id2.kp = id2.kp * id2.taus(aux_i);
%end
%id2.tm = id2.sys.ioDelay;

id2.out = step(id2.sys,time.ext)';
% plot responses
aux_fig = aux_fig + 1;
figure(aux_fig)
plot(time.ext,g123c.out,'r-','LineWidth',2);
hold on
plot(time.ext,sym.out,'b-','LineWidth',2);
hold on
plot(time.ext,stark.out,'g-','LineWidth',2);
hold on
plot(time.ext,jf.out,'m-','LineWidth',2);
hold on
plot(time.ext,ho2.out,'c-','LineWidth',2);
hold on
%plot(time.ext,viteckova2.out,'r:','LineWidth',2);
%hold on
plot(time.ext,id2.out,'y-','LineWidth',2);
hold on
% we need to rescale the original because real step was not unitary
% we could use time.raw and out_raw here if we want...
plot(time.ext,plant.out_ext/max(plant.in),'k-','LineWidth',2);
title('Comparación modelos 2do orden')
%legend(g123c.name,sym.name,stark.name,jf.name,ho2.name,viteckova2.name,id2.name,'Real','Location','SE')
legend(g123c.name,sym.name,stark.name,jf.name,ho2.name,id2.name,'Real','Location','SE')
xlabel('Tiempo (s)')
grid on
saveas(gcf,'comparison_2nd_order','png');

% MSE for 2nd order models
models_2.names = {g123c.name,sym.name,stark.name,jf.name,ho2.name,viteckova2.name,id2.name};
models_2.outs= [g123c.out;sym.out;stark.out;jf.out;ho2.out;viteckova2.out;id2.out];
models_2.number = 7;
models_2.MSE = zeros(models_2.number,1);
%remember to compare against scaled out for unit step
for aux_i = 1:models_2.number
    models_2.MSE(aux_i) = immse(models_2.outs(aux_i,:),plant.out_ext/max(plant.in));
end

[models_2.error,models_2.idx] = min(models_2.MSE);
models_2.best = models_2.names(models_2.idx);

% pay attention to stabilization times for models, even if MSE is low this
% doesnt take into account the whole stabilization time of each system
% should make a longer time and plant out vector to compare better

best.both = [models_1.error,models_2.error];
best.both_names = [models_1.idx ,models_2.idx];
[best.error,best.idx] = min(best.both);

if best.idx == 1
    best.name = models_1.names(best.both_names(best.idx));
end

if best.idx == 2
    best.name = models_2.names(best.both_names(best.idx));
end

best.all = {models_1.MSE, models_1.names';models_2.MSE,models_2.names'};