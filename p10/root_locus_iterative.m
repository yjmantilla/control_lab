clear all;
close all;

%% Modelo 1er Orden
s=tf('s');
kp=1.34;
tau = 78.2;
Gp = kp/(tau*s + 1);
graph_iterations = 1;

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
%un poco a la izquierda(mayor valor absoluto) del polo del sistema
tol = 0.05;%0.0005;
ti = tau + tau*tol;%
comp_zero =  ti*s +1; 
comp_pole = ti*s;

Gc_nok = comp_zero/comp_pole;

bp_limit = 30;
iter = 1;
kcs = [];
bps = [];
sps = [];
tss = [];

tis = [];
chis =[];
wns =[];
Ps=[];
Is=[];
vals = [];
wanteds=[];
reals=[{}];
sp = 0.001:1:10;
ts = 1:1:120;
sps_init = [];
tss_init = [];
for sp = 4%0.001:1:10%30 % subestimated too
   for ts = 99%1:1:120%200 %ts seems to subestimate
        %fprintf('%.2f %.2f\n',sp,ts); 
        chi = sqrt(((log(sp/100))^2)/((pi^2)+(log(sp/100))^2));
        wn = 4/(chi*ts);
        cond5 = wn > 0.0564;
        cond6 = chi > 0.5912;
        if cond5 && cond6
            %fprintf('sol!');

            wanted_root1 = - chi *wn + wn*sqrt(chi*chi-1);
            wanted_root2 = - chi *wn - wn*sqrt(chi*chi-1);
            GpGc_nok = abs(evalfr(Gp*Gc_nok,wanted_root1));
            kc = 1/GpGc_nok;
            %[theta, rho] = cart2pol(real(a), imag(a));
            Gc = kc * Gc_nok;
            sys_r=feedback(Gc*Gp,1);
            real_root = pole(sys_r);
            val = abs(evalfr(Gc*Gp,wanted_root1));
            info = stepinfo(sys_r);
            bp = 100/kc;
            cond1 = info.SettlingTime < 120;
            cond2 = bp >30;
            cond3 = ti >60;
            cond4 = info.Overshoot < 10;
            if (cond1 && cond2 && cond3 && cond4)
            if graph_iterations
                hold on
                step(sys_r);
            end
            fprintf('\n############sol####################!\n')
            %rlocus(Gr)
            
            bps = [bps ; bp];
            tis = [tis ; ti];
            sps = [sps ; info.Overshoot];
            tss = [tss ; info.SettlingTime];
            kcs = [kcs ; kc];
            Ps = [Ps ; kc];
            Is = [Is ; kc/ti];
            chis = [chis;chi];
            wns = [wns;wn];
            vals = [vals;val];
            wanteds = [wanteds;wanted_root1];
            reals = [reals;{real_root}];
            tss_init = [tss_init;ts];
            sps_init = [sps_init;sp];
            end
            
            fprintf('sp %.3f ts %.2f kc %.2f ti %.2f bp %.2f ts %.2f sp %.2f %.3f\n', sp,ts,kc,ti,bp,info.SettlingTime,info.Overshoot,val);

        end
   end
end

[bpmax,bpind] = max(bps);
[spmin,spind] = min(sps);
[tsmin,tsind] = min(tss);
[spmax,spind2] = max(sps);

% Moraleja
% Es dificil poner sp a una planta sin sp (el modelo)
% Pero la planta real real si tiene sp?
% el valor de la parte real deseada se acerca al valor de la parte real
% obtenida
% pero nunca tenemos parte imaginaria?
TPI = table(sps_init,tss_init,bps,sps,tss,chis,wns,vals,kcs,tis,Ps,Is,vals,wanteds,reals);

% analisis conclusiones
% en general ver la diferencia de lo diseñado con lo obtenido en las
% simulaciones


% esto era para ver la trayectoria de la raiz al cambiar la ganancia
% for k_test = 1:0.1:100/bp_limit
%     Gc = k_test * Gc_nok;
%     Gr = feedback(Gp*Gc,1);
%     poles = pole(Gr);
%     p1 = poles(1);
%     p2 = poles(2);
%     re = [re real(poles)];
%     im = [im imag(poles)];
% end
% k_line = 1:0.1:100/bp_limit;

% se elige el 9 {30.5715359503927,0,76.0354860131750,0.715618575909178,0.0564603012892838,1,3.27101654827767,82.1100000000000,3.27101654827767,0.0398370058248407,1,-0.0404040404040404 + 0.0394370275338466i,[-0.0568257128042097;-0.0120126781765073]}