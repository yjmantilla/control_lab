clear all;
close all;
s=tf('s');
%Modelo 1er Orden
kp=1.34;
tau = 78.2;
Gp = kp/(tau*s + 1);
graph_iterations = 0;

sol_re = [];
sol_im = [];
bps = [];
tis = [];
sps = [];
tss = [];
chis = [];
wns = [];
vals =[];
kcs = [];
P1s = [];
I1s = [];

%%% iterativo sin formula
for ti = 60:1 :100
    %disp(ti)
    for kc = 1:0.11:100/30
        bp=100/kc;
        %Gc = kc*(ti*s +1)/ti*s;
        Gc = kc*(1+1/(s*ti));
        Gr = feedback(Gp*Gc,1);
        sol = roots(cell2mat(Gr.den));
        val = evalfr(Gr,sol(1)+i*sol(2));
        r = sol(1)+i*sol(2);
        val=abs((1.34/(78.2*r+1))*((kc*(ti*r+1))/(ti*r)));
        
        wn = sqrt(sol(1)*sol(1) + sol(2)*sol(2));
        chi = -1*sol(1)/wn;
        info = stepinfo(Gr);
        cond1 = info.SettlingTime < 120;
        cond2 = bp > 30;
        cond3 = ti > 60;
        cond4 = info.Overshoot < 10;
        cond5 = wn > 0.0564;
        cond6 = chi > 0.5912;
        if (cond1 && cond2 && cond3 && cond4 && cond5 && cond6)
            if graph_iterations
                hold on
                step(Gr);
            end    
        fprintf('sol:')
        %rlocus(Gr)
        bps = [bps ; bp];
        tis = [tis ; ti];
        sps = [sps ; info.Overshoot];
        tss = [tss ; info.SettlingTime];
        kcs = [kcs ; kc];
        P1s = [P1s ; kc];
        I1s = [I1s ; kc/ti];
        sol_re = [sol_re ; sol(1)];
        sol_im = [sol_im ; sol(2)];
        %phase = sol(2)/sol(1);
        %chi = sqrt(1/(phase*phase+1));
        %wn = -1*sol(1)/chi;
        wns = [wns ; wn];
        chis = [chis; chi];
        vals = [vals; val];
      
        
        end
        %fprintf('%.2f %.2f\n', ti,kc);
        fprintf('kc %.2f ti %.2f bp %.2f ts %.2f sp %.2f\n', kc,ti,bp,info.SettlingTime,info.Overshoot);

        %disp(kc)
    end

end
[bpmax,bpind] = max(bps);
[spmin,spind] = min(sps);
[tsmin,tsind] = min(tss);
[spmax,spind2] = max(sps);

TPI = table(bps,sps,tss,kcs,tis,sol_re,sol_im,chis,wns,vals,P1s,I1s);

% % Seleccionado: (1)
% %% Controlador PI escogido
% sel = 1;
% tp = 78.2;
% wn = wns(sel);
% z = chis(sel);
% kc=(2*z*wn*tp-1)/(kp);
% ti=(kp*kc)/(tp*wn*wn);
% bp=100/kc;
% gc = kc*(1+1/(s*ti));
% g_cl_pi = feedback(Gp*gc,1);
% if graph_iterations
%     title('Respuestas posibles PI')
%     hold off
% end
% figure(2)
% step(g_cl_pi);
% title('Respuesta Escogida PI')
% 
