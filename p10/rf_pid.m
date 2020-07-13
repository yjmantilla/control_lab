clear all
close all
s = tf('s');

%% Modelo 2do orden
kp = 1.34;
Gp_nok = 1/((70.8*s+1)*(22.35*s+1));
Gp = kp*Gp_nok;
alpha = 0.1;

%% Init
td = 5;
factor = 9.1;
mp = 8.5;

%Se realiza un controlador aproximado, para encontrar la fase maxima y la
%frecuencia a la que se da
ti = factor*td;

% PID serie provisional
Gc_xti_divkc = (ti*s+1)*(td*s+1)/(s*(alpha*td*s+1));

% Fase máxima
[mag,phase,w_line] = bode(Gc_xti_divkc);
[Qmax,index] = max(phase);
wmax = w_line(index);

% Margen de fase
z = sqrt(((log(mp/100))^2)/((pi^2)+(log(mp/100))^2));

%Con el valor de z se halla MF tan-1(2z/sqrt(sqrt(1+4*Z^4)-2*z^2)
MF = atan2(2*z,(sqrt(sqrt(1+4*z^4)-2*z^2)))*180/pi;

% Se halla la relacion de diseño con la frecuencia a la que se dio la
% fase max
Rd = 1/(alpha*td*wmax);

% se haya la fase actual;
Fa = -180+ MF- Qmax;

% se realiza el bode a la planta y se mira la frecuencia actual a esa fase
% actual
[mag,phase,w_line] = bode(Gp_nok);
idx = findNearest(phase,Fa);
wactual = w_line(idx);

%se halla los td y ti respectivos del sistema 
tdf = 1/(Rd*wactual*alpha);
tif = factor*tdf;

%se realiza el bode del sistema en red abierta con los nuevos valores de td y ti
%para mirar la fase actual = -180 + MF y a que ganancia esta en el diagrama
%de magnitud

Gc_xti_divkc = (tif*s+1)*(tdf*s+1)/(s*(alpha*tdf*s+1));
G = Gc_xti_divkc*Gp_nok;
[mag,phase,w_line] = bode(G);
idx = findNearest(phase,-180+MF);
K = mag(idx);
w_k = w_line(idx);
K = 1/K;
kc =(K*tif)/(kp);
bp = 100/kc;

% se grafica bode con el fin de mirar que a la fase actual, ya esta en 0dB
Gc_nok = (tif*s+1)*(tdf*s+1)/(tif*s*(alpha*tdf*s+1));
G1 = Gp_nok*Gc_nok*kc*kp;

% se realiza el sistema en red cerrada y se grafica la respeusta al escalon
sys_r = feedback(G1,1);
info = stepinfo(sys_r);
kc_pal = kc*((tif+tdf)/tif);
bp_pal = 100/kc_pal;
ti_pal = tif+tdf;
td_pal = tif*tdf/(tif+tdf);
Gc_pal = kc_pal*(1+(1/(ti_pal*s))+((td_pal*s)/(alpha*td_pal*s+1)));
sys_rpal = feedback(Gp*Gc_pal,1);

info_pal = stepinfo(sys_rpal);
pid.ti = ti_pal;
pid.td = td_pal;
pid.bp = bp_pal;
pid.sp = info_pal.Overshoot;
pid.ts = info_pal.SettlingTime;
sim.Ps = kc_pal;
sim.Is = kc_pal/ti_pal;
sim.Ds = kc_pal*td_pal;
sim.Ns = 1/(alpha*td_pal);
