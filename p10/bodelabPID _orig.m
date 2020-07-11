clear all
close all
s=tf('s');
%Se realiza un controlador aproximado, para encntrar la fase maxima y la
%frecuencia a la que se da
Td=1;
Ti=4*Td;
Gc=(Ti*s+1)*(Td*s+1)/(s*(0.1*Td*s+1));
figure(1)
bode(Gc)
Qmax=50.7;
wmax=3.58;
mp=1;
%z=sqrt(((log(mp/100))^2)/((pi^2)+(log(mp/100))^2));
%Con el valor de chita se halla MF tan-1(2z/sqrt(sqrt(1+4*Z^4)-2*z^2)
z=0.96;
MF=atand(2*z/(sqrt(sqrt(1+4*z^4)-2*z^2)));
% Se halla la relacion de diseño con la frecuencia a la que se dio la
% fase max
Rd=1/(0.1*Td*wmax);
% se haya la fase actual;
Fa=-180+ MF- Qmax;
% se realiza el bode a la planta y se mira la frecuencia actual a esa fase
% actual
kp=1.34;
Gp=1/((70.8*s+1)*(22.35*s+1));
figure(2)
bode(Gp)

wactual=0.134;
%se halla los Td y Ti respectivos del sistema 
Tdf=1/(Rd*wactual*0.1);
Tif=3*Tdf;
% se realiza el bode del sistema en red abierta con los nuevos valores de Td y Ti
%para mirar la fase actual = -180 + MF y a que ganancia esta en el diagrama
%de magnitud
Gc2=(Tif*s+1)*(Tdf*s+1)/(s*(0.1*Tdf*s+1));
G=Gc2*Gp;
figure(3)
bode(G)


KdB=20.7;
% se atenua esa magnitud
K=10^(-KdB/20);
% se calcula la nueva ganancia;
Kc=(K*Tif)/(kp);
% se grafica bode con el fin de mirar que a la fase actual, ya esta en 0dB
Gc2=(Tif*s+1)*(Tdf*s+1)/(Tif*s*(0.1*Tdf*s+1));
G1=Gp*Gc2*Kc*kp;
figure(4)
bode(G1)
% se realiza el sistema en red cerrada y se grafica la respeusta al escalon
sys_r=feedback(G1,1)
figure(5)
step(sys_r)


