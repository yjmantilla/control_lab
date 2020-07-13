clear all;
close all;

%% Modelo 1er Orden
s = tf('s');
kp = 1.34;
tau = 78.2;
Gp = kp/(tau*s + 1);
figure(1)
rlocus(Gp)

%% System poles
plant_poles = pole(Gp);
info = stepinfo(Gp);
orig_ts = info.SettlingTime;
orig_sp = info.Overshoot;

% el root locus es en realidad un conjunto de soluciones
% pero nosotros diseñamos usando su condicion de magnitud
% y forzando una solucion dentro del root locus en particular
% dicha solucion se forma a partir de las relaciones para la
% respuesta caracteristica de segundo orden

%Choose Ti near tau to cancel that pole
tol = 0.05;
ti = tau + tau*tol;
comp_zero =  ti*s +1; 
comp_pole = ti*s;

Gc_nok = comp_zero/comp_pole;

sp = 4;
ts = 99;
chi = sqrt(((log(sp/100))^2)/((pi^2)+(log(sp/100))^2));
wn = 4/(chi*ts);
wanted_root1 = - chi *wn + wn*sqrt(chi*chi-1);
wanted_root2 = - chi *wn - wn*sqrt(chi*chi-1);
GpGc_nok = abs(evalfr(Gp*Gc_nok,wanted_root1));
kc = 1/GpGc_nok;
Gc = kc * Gc_nok;
sys_r=feedback(Gc*Gp,1);
figure(2)
rlocus(sys_r);
real_root = pole(sys_r);
val = abs(evalfr(Gc*Gp,wanted_root1));
info = stepinfo(sys_r);
pi.bp = 100/kc;
pi.kc = kc;
pi.ti = ti;
pi.sp = info.Overshoot;
pi.ts = info.SettlingTime;
sim.P = kc;
sim.I = kc/ti;