close all;
clear all;
s=tf('s');
graph_iterations = 1; % 0 faster, 1 graphs posible solutions

%% Modelo 1er orden
kp=1.34; % Ganancia de la planta
tm=8.4; % Tiempo Muerto
tp=78.2; % Constante de tiempo de la planta
g1= (kp*exp(-tm*s))/(tp*s+1);
g1_notm =(kp)/(tp*s+1); % Funcion de la planta sin tiempo muerto

%% Diseño de controlador PI

bps = [];
tis = [];
mps = [];
tss = [];
wns = [];
kcs = [];
P1s = [];
I1s = [];

% Metodo de binomio estandarizado
z=1;
if graph_iterations
    figure(1)
    hold on
end

for wn=0.01:0.0005:0.03
kc=(2*z*wn*tp-1)/(kp);
ti=(kp*kc)/(tp*wn*wn);
bp=100/kc;
gc = kc*(1+1/(s*ti));
g_cl_pi = feedback(g1_notm*gc,1);

info = stepinfo(g_cl_pi);
cond1 = info.SettlingTime < 120;
cond2 = bp > 30;
cond3 = ti > 60;
cond4 = info.Overshoot < 10;

if (cond1 && cond2 && cond3 && cond4)
    if graph_iterations
        step(g_cl_pi);
    end    
    fprintf('sol:')
    wns = [wns ; wn];
    bps = [bps ; bp];
    tis = [tis ; ti];
    mps = [mps ; info.Overshoot];
    tss = [tss ; info.SettlingTime];
    kcs = [kcs ; kc];
    P1s = [P1s ; kc];
    I1s = [I1s ; kc/ti];
end
disp(wn)
end
TPI = table(wns,bps,mps,tss,kcs,tis,P1s,I1s);

%% Controlador PI escogido
wn = 0.0235;
z = 1;
kc=(2*z*wn*tp-1)/(kp);
ti=(kp*kc)/(tp*wn*wn);
bp=100/kc;
gc = kc*(1+1/(s*ti));
g_cl_pi = feedback(g1_notm*gc,1);
if graph_iterations
    title('Respuestas posibles PI')
    hold off
end
figure(2)
step(g_cl_pi);
title('Respuesta Escogida PI')

%% Modelo 2do orden

kp=1.34; % ganancia de la planta
g22 = kp/((70.8*s+1)*(22.35*s+1));

% Las siguientes constantes resulta de volver el denominador de g22 monico
k = 8.468e-4; 
a = 0.0588;
b = 6.319e-4;
g2=tf([k],[1 a b]);

%% PID
b_ps = [];
t_is = [];
m_ps = [];
t_ss = [];
w_os = [];
k_cs = [];
t_ds = [];
P2_s = [];
I2_s = [];
D2_s = [];
N2_s = [];
alpha = 0.1; % Relacionado al retardo de la accion derivativa
% Metodo binomial estandarizado 1 3 3 1

if graph_iterations
    figure(3)
    hold on
end

for wo=0.0005:0.0005:0.05;
    k_c=(3*wo*wo-b)/(k); %2.15 ITAE
    t_i=(k_c*k)/(wo*wo*wo);
    t_d=((3*wo)-a)/(k*k_c); %1.75 ITAE
    b_p=100/k_c;
    c = k_c*(1+(1/(t_i*s))+((t_d*s)/(alpha*t_d*s+1)));
    %c = k_c*(1+(1/(t_i*s))+((t_d*s))); % Si se quiere sin el retardo
    gp_2=tf(1.34,[1582.38 93.15 1]);
    sys=c*gp_2;
    sys_r=feedback(sys,1);
    info = stepinfo(sys_r);
    cond1 = info.SettlingTime < 120;
    cond2 = b_p > 30;
    cond3 = t_i > 60;
    cond4 = info.Overshoot < 10;

    if (cond1 && cond2 && cond3 && cond4)
        if graph_iterations
            step(sys_r);
        end
        
        fprintf('sol:')
        w_os = [w_os ; wo];
        b_ps = [b_ps ; b_p];
        t_is = [t_is ; t_i];
        m_ps = [m_ps ; info.Overshoot];
        t_ss = [t_ss ; info.SettlingTime];
        k_cs = [k_cs ; k_c];
        t_ds = [t_ds ; t_d];
        P2_s = [P2_s ; k_c];
        I2_s = [I2_s ; k_c/t_i];
        D2_s = [D2_s ; k_c*t_d];
        N2_s = [N2_s ; 1/(alpha*t_d)];
    end

    disp(wo)

end

TPID = table(w_os,b_ps,m_ps,t_ss,k_cs,t_is,t_ds,P2_s,I2_s,D2_s,N2_s);

%% Controlador PID escogido
wo = 0.027;
k_c=(3*wo*wo-b)/(k); %2.15 ITAE
t_i=(k_c*k)/(wo*wo*wo);
t_d=((3*wo)-a)/(k*k_c); %1.75 ITAE
b_p=100/k_c;
c = k_c*(1+(1/(t_i*s))+((t_d*s)/(alpha*t_d*s+1)));
%c = k_c*(1+(1/(t_i*s))+((t_d*s))); % Si se quiere sin el retardo
gp_2=tf(1.34,[1582.38 93.15 1]);
sys=c*gp_2;
sys_r=feedback(sys,1);
if graph_iterations
    title('Respuestas posibles PID')
    hold off
end
figure(4)
step(sys_r);
title('Respuesta Escogida PID')
