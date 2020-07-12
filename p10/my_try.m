clear all
close all
s=tf('s');
%% 2nd order model
plant.kp = 1.34;
plant.taus = [70.8;22.35];
plant.G_nok = 1/((plant.taus(1)*s+1)*(plant.taus(2)*s+1));
plant.G = plant.kp*plant.G_nok;
alpha = 0.1;

%Input Parameters
param.td = 1;
param.factor = 10;
param.sp = 5;

% Approximate Controller
param.ti = param.factor*param.td;

% You can check that the non complex parts of Gc don't influence the phase
% kc = 2.3;
% Gc_approx = kc*(ti*s+1)*(td*s+1)/(s*ti*(alpha*td*s+1));
control.approx.G_xti_divk = (param.ti*s+1)*(param.td*s+1)/(s*(alpha*param.td*s+1));
[bode.mag,bode.phase,bode.w]= bode(control.approx.G_xti_divk);
[phase.max.phase,index]= max(bode.phase);
phase.max.w = bode.w(index);
clear bode index

% Find a z given sp to input in the phase margin formula
z=sqrt(((log(param.sp/100))^2)/((pi^2)+(log(param.sp/100))^2));
phase.margin = atan2(2*z,(sqrt(sqrt(1+4*z^4)-2*z^2)))*180/pi;

% Design Relationship 
param.dr = (1/(alpha*param.td))/phase.max.w;

% The controller changes the phase margin of the system by adding its phase
% phase.margin = phase.unity + phase.max(by controller) + 180
% phase.unity is the unity-gain phase (or gain crossover phase)
phase.unity.phase = phase.margin - 180 - phase.max.phase;

% Check the frequency of the plant at unity-gain 
% its gain crossover frequency
[bode.mag,bode.phase,bode.w]= bode(plant.G_nok);
idx = findNearest(bode.phase,phase.unity.phase);
phase.unity.w = bode.w(idx);
clear bode idx

% Design the real td using the previous design relationship
control.series.td = 1/(param.dr*phase.unity.w*alpha);
control.series.ti = param.factor*control.series.td;

% Check designed controller + plant at the desired unity-gain phase
% phase.margin = closed_sys.phase.unity + 180 
control.G_xti_divk =(control.series.ti*s+1)*(control.series.td*s+1)/(s*(alpha*control.series.td*s+1));
closed.G_nok = control.G_xti_divk*plant.G_nok;
[mag,phase,w_line] = bode(closed.G_nok);
idx = findNearest(phase,-180+phase.margin);
K = mag(idx);
K = 1/K;
kc=(K*tif)/(plant.kp);

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