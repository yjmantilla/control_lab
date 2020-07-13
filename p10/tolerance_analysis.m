% %% Se necesita por lo menos matlab 2014b
% para aplicar la funcion round(x,n) (de 2 parametros)
close all;
clear all;
s=tf('s');
graph_iterations = 0; % 0 faster, 1 graphs posible solutions
tol_step = 0.1;
req.ts = 120;
req.ov = 10;
req.ti = 60;
req.bp = 30;
n = 3;

%% Plant Models
kp=1.34;
Gp2 = kp/((70.8*s+1)*(22.35*s+1));
Gp1 = kp/(78.2*s + 1);

alpha = 0.1;
%% Designed PIDs

% P9
p9.pi.kc = 1.997;
p9.pi.ti = 61.95;
p9.pi.td = 0;

p9.pid.kc = 1.84;
p9.pid.ti = 79.01;
p9.pid.td = 14.28;

%P10
p10.pi.kc = 3.7;
p10.pi.ti = 82.11;
p10.pi.td = 0;

p10.pid.kc = 3.07;
p10.pid.ti = 71.7;
p10.pid.td = 6.4;

%All
All = [p9.pi.kc , p9.pi.ti  , p9.pi.td;
p9.pid.kc , p9.pid.ti , p9.pid.td;
p10.pi.kc , p10.pi.ti , p10.pi.td;
p10.pid.kc, p10.pid.ti, p10.pid.td];

source = {'AP_PI';'AP_PID';'LR_PI';'RF_PID'};
s_nums = [1,2,3,4];
% You may analize how the controllers perform on differents orders of the
% system with the following array ie [1,1,1,1] , [2,2,2,2]
% or as originally designed [1,2,1,2] 
orders = [1,1,1,1];

%% Setup Parameters
kcs = [];
tis = [];
tds = [];
bps = [];
tss = [];
ovs = [];
Ps = [];
Is = [];
Ds = [];
Ns = [];
sources = [];
plant_orders =[];

for i =1:length(All)
kc = All(i,1);
ti = All(i,2);
td = All(i,3);

kc_vec = linspace(kc - kc*0.15,kc + kc*0.15,n);
ti_vec = linspace(ti - ti*0.15,ti + ti*0.15,n);
td_vec = linspace(td - td*0.15,td + td*0.15,n);

kc_vec = [kc_vec, kc];
ti_vec = [ti_vec, ti];
td_vec = [td_vec, td];


%% Iteration
for kc = kc_vec
    for ti = ti_vec
    for td = td_vec
    if orders(i) == 2
        Gp = Gp2;
    end
    if orders(i) == 1
        Gp = Gp1;
    end
    fprintf('kc %.2f ti %.2f td %2.f\n',kc,ti,td);
    bp = 100/kc;
    Gc = kc*(1+1/(s*ti)+((td*s)/(alpha*td*s+1)));
    Gr = feedback(Gp*Gc,1);
    monic = Gr.den{1}(1);
    Gr_monic.num = Gr.num{1}/monic;
    Gr_monic.den = Gr.den{1}/monic;
    Gr_monic.num = tf(Gr.num{1}/monic,1);
    Gr_monic.den = tf(Gr.den{1}/monic,1);
    info = stepinfo(Gr);
    ts = info.SettlingTime;
    ov = info.Overshoot;
    cond = ts <= req.ts && bp >= req.bp && ov <= req.ov && ti >= req.ti;
    if cond
    kcs = [kcs ; kc];
    tis = [tis ; ti];
    tds = [tds ; td];
    bps = [bps ; bp];
    tss = [tss ; ts];
    ovs = [ovs ; ov];
    Ps = [Ps ; kc];
    Is = [Is ; kc/ti];
    Ds = [Ds ; kc*td];
    plant_orders =[plant_orders;orders(i)];
    sources = [sources; s_nums(i)];
    if td ~= 0
        Ns = [Ns ; 1/(alpha*td)];
    else
        Ns = [Ns ; 0];
    end
    if graph_iterations
        figure(1)
        hold on
        step(Gr);
    end
    
    end
    end
    end
end
end
% if matlab is less than R2014b comment the following rounds
kcs = round(kcs,2);
tis = round(tis,2);
tds = round(tds,2);
bps = round(bps,2);
tss = round(tss,2);
ovs = round(ovs,2);
Ps  = round(Ps ,2);
Is  = round(Is ,2);
Ds  = round(Ds ,2);
Ns  = round(Ns ,2);
TPID = table(sources,plant_orders,kcs,tis,tds,bps,tss,ovs,Ps,Is,Ds,Ns);
TPID = unique(TPID);
TPID = sortrows(TPID,'tss','ascend');