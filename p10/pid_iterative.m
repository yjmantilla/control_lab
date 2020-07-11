clear all;
close all;
s=tf('s');
kp=1.34;
Gp = kp/(78.2*s + 1);
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
Ps = [];
Is = [];
Ds = [];
Ns = [];
tds = [];
alpha = 0.1;
%%% iterativo sin formula
for ti = 60  :1 :100
    %disp(ti)
    for kc = 1:0.1:100/30
            for td = ti/20:1:ti/5

            bp=100/kc;
            %Gc = kc*(ti*s +1)/ti*s;
            Gc = kc*(1+(1/(ti*s))+((td*s)/(alpha*td*s+1)));
            Gr = feedback(Gp*Gc,1);

            info = stepinfo(Gr);
            cond1 = info.SettlingTime < 120;
            cond2 = bp > 30;
            cond3 = ti > 60;
            cond4 = info.Overshoot < 10;
            cond5 = 1;%wn > 0.0564;
            cond6 = 1;%chi > 0.5912;
                if (cond1 && cond2 && cond3 && cond4 && cond5 && cond6)
                            if graph_iterations
                                hold on
                                step(Gr);
                            end    
                    fprintf('sol:')
                    %rlocus(Gr)
                    bps = [bps ; bp];
                    tis = [tis ; ti];
                    tds = [tds ; td];
                    sps = [sps ; info.Overshoot];
                    tss = [tss ; info.SettlingTime];
                    kcs = [kcs ; kc];
                    Ps = [Ps ; kc];
                    Is = [Is ; kc/ti];
                    Ds = [Ds ; kc*td];
                    Ns = [Ns ; 1/(alpha*td)];
                  

                end
            fprintf('kc %.2f ti %.2f td %.2f bp %.2f ts %.2f sp %.2f\n', kc,ti,td,bp,info.SettlingTime,info.Overshoot);
            end

        %disp(kc)
    end

end
[bpmax,bpind] = max(bps);
[spmin,spind] = min(sps);
[tsmin,tsind] = min(tss);
[spmax,spind2] = max(sps);

TPID = table(bps,sps,tss,kcs,tis,sol_re,sol_im,chis,wns,vals,Ps,Is,Ds,Ns);

% Seleccionado: (1)
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
