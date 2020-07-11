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
kcs = [];
tis = [];
tds = [];

bps = [];
sps = [];
tss = [];

Ps = [];
Is = [];
Ds = [];
Ns = [];
j=1;
for mp=0.001:0.5:10
   for td=1:10
       for factor=60/td:1:20
        fprintf('\n%.3f %.2f %.2f ', mp,td,factor);
        %Se realiza un controlador aproximado, para encntrar la fase maxima y la
        %frecuencia a la que se da
        ti=factor*td;
        Gc=(ti*s+1)*(td*s+1)/(s*(alpha*td*s+1));
        [mag,phase,w_line]= bode(Gc);
        [Qmax,index]= max(phase);
        wmax=w_line(index);
        %Qmax=50.7;
        %wmax=3.58;
        %0.005;
        z=sqrt(((log(mp/100))^2)/((pi^2)+(log(mp/100))^2));
        %Con el valor de chita se halla MF tan-1(2z/sqrt(sqrt(1+4*Z^4)-2*z^2)
        %z=0.96;
        MF=atan2(2*z,(sqrt(sqrt(1+4*z^4)-2*z^2)))*180/pi;
        % Se halla la relacion de diseño con la frecuencia a la que se dio la
        % fase max
        Rd=1/(0.1*td*wmax);
        % se haya la fase actual;
        Fa=-180+ MF- Qmax;
        % se realiza el bode a la planta y se mira la frecuencia actual a esa fase
        % actual
        [mag,phase,w_line]= bode(Gp_nok);
        idx = findNearest(phase,Fa);
        wactual = w_line(idx);
        %wactual=0.134;
        %se halla los td y ti respectivos del sistema 
        tdf=1/(Rd*wactual*0.1);
        tif=factor*tdf;%3?
        cond3 = tif > 60;
        if not(cond3)
            continue
        end
        % se realiza el bode del sistema en red abierta con los nuevos valores de td y ti
        %para mirar la fase actual = -180 + MF y a que ganancia esta en el diagrama
        %de magnitud
        Gc2=(tif*s+1)*(tdf*s+1)/(s*(alpha*tdf*s+1));
        G=Gc2*Gp_nok;
        [mag,phase,w_line] = bode(G);

        idx = findNearest(phase,Fa);
        K = mag(idx);
        %KdB=20.7;
        % se atenua esa magnitud
        %K=10^(-KdB/20);
        % se calcula la nueva ganancia;
        %kc=(tif)/(K*kp);
        kc=(K*tif)/(kp);
        bp= 100/kc;
        cond2 = bp > 30;
        if not (cond2)
            continue
        end
        % se grafica bode con el fin de mirar que a la fase actual, ya esta en 0dB
        Gc2=(tif*s+1)*(tdf*s+1)/(tif*s*(alpha*tdf*s+1));
        G1=Gp_nok*Gc2*kc*kp;

        %bode(G1)
        % se realiza el sistema en red cerrada y se grafica la respeusta al escalon
        sys_r=feedback(G1,1);
        info = stepinfo(sys_r);
        
        cond1 = info.SettlingTime < 120;
        
        
        cond4 = info.Overshoot < 10;
        %step(sys_r)
        if (cond1 && cond2 && cond3 && cond4)
            if graph_iterations
                hold on
                step(sys_r);
            end
            fprintf('sol!!!!!!!!!!!!!')
            %rlocus(Gr)
            
            bps = [bps ; bp];
            tis = [tis ; tif];
            tds = [tds ; tdf];
            sps = [sps ; info.Overshoot];
            tss = [tss ; info.SettlingTime];
            kcs = [kcs ; kc];
            Ps = [Ps ; kc];
            Is = [Is ; kc/tif];
            Ds = [Ds ; kc*tdf];
            Ns = [Ns ; 1/(alpha*td)];
            
        end
       %fprintf('%.3f %.2f %.2f kc %.2f ti %.2f td %.2f bp %.2f ts %.2f sp %.2f\n', mp,td,factor,kc,ti,td,bp,info.SettlingTime,info.Overshoot);
       
       end
   
   end
end

[bpmax,bpind] = max(bps);
[spmin,spind] = min(sps);
[tsmin,tsind] = min(tss);
[spmax,spind2] = max(sps);

TPID = table(bps,sps,tss,kcs,tis,tds,Ps,Is,Ds,Ns);

for i=1:length(bps)
    %i=1;%select individual graph
    Gc = kcs(i)*(1+(1/(tis(i)*s))+((tds(i)*s)/(alpha*tds(i)*s+1)));
    %bode(G1)
    % se realiza el sistema en red cerrada y se grafica la respeusta al escalon
    sys_r=feedback(Gp*Gc,1);
    hold on
    step(sys_r);
end

% analisis conclusiones
% en general ver la diferencia de lo diseñado con lo obtenido en las
% simulaciones
