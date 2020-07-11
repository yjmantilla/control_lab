clear all
close all
s=tf('s');

%Se realiza un controlador aproximado, para encontrar la fase maxima y la
%frecuencia a la que se da

Td=1;
Ti=4*Td;
kp=1.34;
Gc=(Ti*s+1)*(Td*s+1)/(s*(0.1*Td*s+1));

[mag,phase,w_line] = bode(Gc);
phase_max = struct;
[phase_max.val,phase_max.index] = max(phase);
phase_max.w = w_line(phase_max.index);

mp=1;
z=sqrt(((log(mp/100))^2)/((pi^2)+(log(mp/100))^2));

%Con el valor de chita se halla MF tan-1(2z/sqrt(sqrt(1+4*Z^4)-2*z^2)
MF=atan2(2*z,(sqrt(sqrt(1+4*z^4)-2*z^2)))*180/pi;

% Se halla la relacion de diseño con la frecuencia a la que se dio la
% fase max

Rd=1/(0.1*Td*phase_max.w);

% se haya la fase actual;
Fa=-180+ MF- phase_max.val;

% se realiza el bode a la planta y se mira la frecuencia actual a esa fase
% actual
Gp.k = 1.34;
Gp.no_k=1/((70.8*s+1)*(22.35*s+1));
[mag,phase,w_line] = bode(Gp.no_k);

f_line = w_line*180/pi;
idx = findNearest(phase,Fa);
wactual = w_line(idx);

%se halla los Td y Ti respectivos del sistema 
Tdf=1/(Rd*wactual*0.1);
Tif=4*Tdf;

% se realiza el bode del sistema en red abierta con los nuevos valores de Td y Ti
%para mirar la fase actual = -180 - MF y a que ganancia esta en el diagrama
%de magnitud
Gc2=(Tif*s+1)*(Tdf*s+1)/(s*(0.1*Tdf*s+1));
G=Gc2*Gp.no_k;
%figure(3)
[mag,phase,w_line] = bode(G);
f_line = w_line*180/pi;
idx = findNearest(phase,Fa);
K = mag(idx);
KdB = 20*log10(K);
Kc=Tif/(K*1.34);
% se grafica bode con el fin de mirar que a la fase actual, ya esta en 0dB
G1=Gc2*Gp.no_k*Kc*Gp.k;
figure(4)
bode(G1)
% se realiza el sistema en red cerrada y se grafica la respeusta al escalon
sys_r=feedback(G1,1);
figure(5)
step(sys_r)


