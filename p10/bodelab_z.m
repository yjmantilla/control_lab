clear all
close all
s=tf('s');
%% Modelo 2do orden
kp=1.34;
Gp_nok=1/((70.8*s+1)*(22.35*s+1));
Gp=kp*Gp_nok;
alpha = 0.1;
%%Init
graph_iterations = 0;
tis_init=[];
tds_init=[];
mps_init = [];
kcs = [];
tis = [];
tds = [];
kcs_pal = [];
tds_pal =[];
tis_pal= [];
bps = [];
sps = [];
tss = [];
bps_pal = [];
sps_pal = [];
tss_pal = [];
Ps = [];
Is = [];
Ds = [];
Ns = [];
j=1;
factors=[];
sol_bool = 0;
lim = 60;
mp =0;
td=7;
factor = 8;
mp = 3.1;
for z=0.6:0.05:1%0.001:1:10
   for td=1:1:20%7%6.9:0.05:7.1%7%1:1:20
       for factor=50%lim/td:0.5:30%8%7.9:0.1:8.1%lim/td:0.5:30
        mp=0;
        fprintf('\n%.3f %.2f %.2f ', z,td,factor);
        %Controlador aproximado PID serie
        ti=factor*td;
        Gc_xti_divkc=(ti*s+1)*(td*s+1)/(s*(alpha*td*s+1));%ti dividiendo?
        
        % Max phase by this controller
        [mag,phase,w_line]= bode(Gc_xti_divkc);
        [phase_max,index]= max(phase);
        wmax = w_line(index);
        
        % Look the desired phase margin using z
        MF=atan2(2*z,(sqrt(sqrt(1+4*z^4)-2*z^2)))*180/pi;
        % Se halla la relacion de diseño con la frecuencia a la que se dio la
        % fase max
        % Design relation = alpha td pole/wmax
        Rd=(1/(alpha*td))/wmax;
        % se haya la fase actual, el controlador agrega a la fase actual
        %ya que se multiplica por la planta
        % Fa + phase_max + 180 = MF
        Fa=-180+ MF- phase_max;
        % se realiza el bode a la planta y se mira la frecuencia actual a esa fase
        % actual
        [mag,phase,w_line]= bode(Gp_nok);
        idx = findNearest(phase,Fa);
        wactual = w_line(idx);%0.37;%0.231
        %wactual=0.134;
        %se halla los td y ti respectivos del sistema 
        tdf=1/(Rd*wactual*alpha);
        tif=factor*tdf;%3?
        cond3 = tif > 60;
        if not(cond3)
            continue
        end
        % se realiza el bode del sistema en red abierta con los nuevos valores de td y ti
        %para mirar la fase actual = -180 + MF y a que ganancia esta en el diagrama
        %de magnitud
        Gc_xti_divkc =(tif*s+1)*(tdf*s+1)/(s*(alpha*tdf*s+1));
        G=Gc_xti_divkc*Gp_nok;
        [mag,phase,w_line] = bode(G);

        idx = findNearest(phase,-180+MF);
        K = mag(idx);
        K = 1/K;
        %KdB = -1*20*log10(K);
        %KdB=20.7;
        % se atenua esa magnitud
        %K=10^(KdB/20);
        % se calcula la nueva ganancia;
        %kc=(tif)/(K*kp);
        kc=(K*tif)/(kp);
        bp= 100/kc;
        cond2 = bp > 30;
        %if not (cond2)
            %continue
        %end
        % se grafica bode con el fin de mirar que a la fase actual, ya esta en 0dB
        Gc_nok=(tif*s+1)*(tdf*s+1)/(tif*s*(alpha*tdf*s+1));
        G1=Gp_nok*Gc_nok*kc*kp;

        %bode(G1)
        % se realiza el sistema en red cerrada y se grafica la respeusta al escalon
        sys_r=feedback(G1,1);
        info = stepinfo(sys_r);
        kc_pal = kc*((tif+tdf)/tif);
        bp_pal = 100/kc_pal;
        ti_pal = tif+tdf;
        td_pal = tif*tdf/(tif+tdf);

        cond1 = info.SettlingTime < 120;
        cond5 = bp_pal > 30;
        cond6 = ti_pal > 60;
        
        cond4 = info.Overshoot < 10;
        %step(sys_r)
        if (cond1 && cond2 && cond3 && cond4 && cond5 && cond6)
            if graph_iterations
                hold on
                step(sys_r);
            end
            fprintf('sol!!!!!!!!!!!!!')
            %rlocus(Gr)
        end
            %Gc_pal = kc_pal*(1+(1/(tis_pal*s))+((tds_pal*s)/(alpha*tds_pal*s+1)));
            %sys_r=feedback(Gp*Gc_pal,1);
            %hold on
            %info_pal = stepinfo(sys_r);
            %info_pal.Overshoot;
            %info_pal.SettlingTime;
            tis_init = [tis_init;ti];
            factors= [factors;factor];
            tds_init = [tds_init;td];
            mps_init = [mps_init;mp];
            bps = [bps ; bp];
            tis = [tis ; tif];
            tds = [tds ; tdf];
            sps = [sps ; info.Overshoot];
            tss = [tss ; info.SettlingTime];
            kcs = [kcs ; kc];
            kc_pal = kc*((tif+tdf)/tif);
            ti_pal = tif+tdf;
            td_pal = tif*tdf/(tif+tdf);
            kcs_pal = [kcs_pal;kc_pal];
            tis_pal = [tis_pal;ti_pal];
            tds_pal = [tds_pal;td_pal];
            
            Ps = [Ps ; kc_pal];
            Is = [Is ; kc_pal/ti_pal];
            Ds = [Ds ; kc_pal*td_pal];
            Ns = [Ns ; 1/(alpha*td_pal)];
            %bps_pal = [bps_pal ; bp_pal];
            %sps_pal = [sps_pal ; info_pal.Overshoot];
            %tss_pal = [tss_pal ;info_pal.SettlingTime];
            sol_bool = 0;
            %break
            
        end
       %fprintf('%.3f %.2f %.2f kc %.2f ti %.2f td %.2f bp %.2f ts %.2f sp %.2f\n', mp,td,factor,kc,ti,td,bp,info.SettlingTime,info.Overshoot);
       if sol_bool
              break
       end
       %end
          if sol_bool
              break
          end
   end
          if sol_bool
              break
          end
end

[bpmax,bpind] = max(bps);
[spmin,spind] = min(sps);
[tsmin,tsind] = min(tss);
[spmax,spind2] = max(sps);
bps_pal = 100./kcs_pal;
%bps_pal = bps_pal';
TPID = table(bps,sps,tss,kcs,tis,tds,Ps,Is,Ds,Ns,kcs_pal,tis_pal,tds_pal,bps_pal,mps_init,tds_init,factors,tis_init);
% 
% for i=1:length(bps)
%     %i=1;%select individual graph
%     Gc = kcs_pal(i)*(1+(1/(tis_pal(i)*s))+((tds_pal(i)*s)/(alpha*tds_pal(i)*s+1)));
%     %bode(G1)
%     % se realiza el sistema en red cerrada y se grafica la respeusta al escalon
%     sys_r=feedback(Gp*Gc,1);
%     hold on
%     step(sys_r);
% end

% analisis conclusiones
% en general ver la diferencia de lo diseñado con lo obtenido en las
% simulaciones
% MFS = [];
% ZS = [];
% for z = 0.1:0.01:1
%     
%     MF=atan2(2*z,(sqrt(sqrt(1+4*z^4)-2*z^2)))*180/pi;
%     MFS=[MFS;MF];
%     ZS = [ZS;z];
%     disp(MF)
% end