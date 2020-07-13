clear all
close all
s=tf('s');
%% Modelo 2do orden
kp=1.34;
Gp_nok=1/((70.8*s+1)*(22.35*s+1));
Gp=kp*Gp_nok;
alpha = 0.1;
%%Init
graph_iterations = 1;
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
lim_ti = 60;
lim_bp = 30;
lim_ts = 123;
lim_ov = 12;
mp =0;
td=7;
factor = 8;
mp = 3.1;
for mp=8.5%:0.01:8.6%8:0.5:9.5%0.001:1:10 %3.1%:0.05:3.2%
   for td=5%1:1:20%7%6.9:0.05:7.1%7%
       for factor=9.1%9.135%:0.0000002:9.1361%lim_ti/td:0.5:15%8%7.9:0.1:8.1%
        fprintf('\n%.3f %.2f %.7f ', mp,td,factor);
        %Se realiza un controlador aproximado, para encntrar la fase maxima y la
        %frecuencia a la que se da
        ti=factor*td;
        %pid serie
        Gc_xti_divkc=(ti*s+1)*(td*s+1)/(s*(alpha*td*s+1));%ti dividiendo?
        [mag,phase,w_line]= bode(Gc_xti_divkc);
        [Qmax,index]= max(phase);
        wmax=w_line(index);
        Qmax=52.98;
        wmax=0.6776;
        %0.005;
        z=sqrt(((log(mp/100))^2)/((pi^2)+(log(mp/100))^2));
        %z=0;
        %Con el valor de chita se halla MF tan-1(2z/sqrt(sqrt(1+4*Z^4)-2*z^2)
        %z=0.96;
        %z=0;
        MF=atan2(2*z,(sqrt(sqrt(1+4*z^4)-2*z^2)))*180/pi;
        % Se halla la relacion de diseño con la frecuencia a la que se dio la
        % fase max
        Rd=1/(alpha*td*wmax);
        % se haya la fase actual;
        
        Fa=-180+ MF- Qmax;
        % se realiza el bode a la planta y se mira la frecuencia actual a esa fase
        % actual
        [mag,phase,w_line]= bode(Gp_nok);
        idx = findNearest(phase,Fa);
        wactual = w_line(idx);%0.37;%0.231
        %wactual = 0.476;
        %wactual=0.134;
        %se halla los td y ti respectivos del sistema 
        tdf=1/(Rd*wactual*alpha);
        tif=factor*tdf;%3?
        tdf = 7.1;
        tif = 64.6;%factor*tdf;
        cond3 = tif > lim_ti;
        if not(cond3)
            continue
        end
        % se realiza el bode del sistema en red abierta con los nuevos valores de td y ti
        %para mirar la fase actual = -180 + MF y a que ganancia esta en el diagrama
        %de magnitud
        Gc_xti_divkc =(tif*s+1)*(tdf*s+1)/(s*(alpha*tdf*s+1));
        G=Gc_xti_divkc*Gp_nok;
        [mag,phase,w_line] = bode(G);
        %dcm = datacursormode(gcf) ;
        idx = findNearest(phase,-180+MF);
        K = mag(idx);
        w_k = w_line(idx);
        K = 1/K;
        KdB = 24.83;
        %KdB=20*log10(mag(idx));%.392;
        % se atenua esa magnitud
        K2=10^(-KdB/20);
        %tdf=round(tdf,1);
        %tif=round(tif,1);
        % se calcula la nueva ganancia;
        %kc=(tif)/(K*kp);
        %kc=(K*tif)/(kp);
        kc=(K2*tif)/(kp);
        kc = 2.77;
        bp= 100/kc;
        cond2 = bp > lim_bp;
        if not (cond2)
            continue
        end

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
        Gc_pal = kc_pal*(1+(1/(ti_pal*s))+((td_pal*s)/(alpha*td_pal*s+1)));
        sys_rpal=feedback(Gp*Gc_pal,1);

        cond1 = info.SettlingTime < lim_ts;
        cond5 = bp_pal > lim_bp;
        cond6 = ti_pal > lim_ti;
        
        cond4 = info.Overshoot < lim_ov;
        %step(sys_r)
            if (cond1 && cond2 && cond3 && cond4 && cond5 && cond6)
            if graph_iterations
                %figure(1)
                %hold on
                %step(sys_r);
                %figure(2)
                %hold on
                step(sys_rpal);
            end
            fprintf('sol!!!!!!!!!!!!!')
            %rlocus(Gr)
        %end
            %hold on
            info_pal = stepinfo(sys_rpal);
            %sp_pal = info_pal.Overshoot;
            %sp_ts = info_pal.SettlingTime;
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

            kcs_pal = [kcs_pal;kc_pal];
            tis_pal = [tis_pal;ti_pal];
            tds_pal = [tds_pal;td_pal];
            
            Ps = [Ps ; kc_pal];
            Is = [Is ; kc_pal/ti_pal];
            Ds = [Ds ; kc_pal*td_pal];
            Ns = [Ns ; 1/(alpha*td_pal)];
            bps_pal = [bps_pal ; bp_pal];
            sps_pal = [sps_pal ; info_pal.Overshoot];
            tss_pal = [tss_pal ;info_pal.SettlingTime];
            sol_bool = 0;
            %break
            end
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
%bps_pal = 100./kcs_pal;
%bps_pal = bps_pal';
TPID = table(bps,sps,tss,kcs,tis,tds,Ps,Is,Ds,Ns,kcs_pal,tis_pal,tds_pal,bps_pal,sps_pal,tss_pal,mps_init,tds_init,factors,tis_init);
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