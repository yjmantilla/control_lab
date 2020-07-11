k=12.5; %gain
ng=1; %num of G
dg=poly([0 -2.5 -5]); %den of G
% create a logarithmically spaced vector 1e-1=0.1 to 1e1=10 with 100 points
w=logspace(-1,1,100)';

% Use any of the following 
bode(k*ng,dg,w);grid % From poly vector
sys=tf(k*ng,dg);bode(sys,w); grid ;% From dynamic system object
bode(tf(k*ng,dg),w); grid ;% From transfer function

% Or use symbolic s to make the transfer function
s=tf('s'); g = 12.5 /s/(s+2.5)/(s+5); bode(g,w); grid;
[mag,phase] = bode(g,w);
mag = 20*log10(squeeze(mag));
phase = squeeze(phase);
idx_cg = findNearest(0,mag); % cruce de ganancia
idx_aba = findNearest(-3,mag); %ancho de banda en red abierta
idx_cf = findNearest(-180,phase);
w_cg = w(idx_cg);
wm_aba = w(idx_aba);
wm_cf = w(idx_cf);

figure()

%[Gm,Pm,Wcg,Wcp] = margin(SYS)
[Gm,Pm,Wcg,Wcp] = margin(g);


sysab=tf(k*ng,dg); % define open loop system
sysce=feedback(sysab,1); % apply negative unitary feedback 
figure(); bode(sysce,w); grid  % plot FR of closed loop system
[magF,phaseF]= bode(sysce,w);
magF = 20*log10(squeeze(magF));
phaseF = squeeze(phaseF);
Mr = max(magF);
idx_rF = findNearest(0,Mr); % resonant index
idx_bwF = findNearest(-3,magF); %ancho de banda en red cerrada
wr = w(idx_rF);
wabc = w(idx_bwF);

% low pass quality of feedback control systems

t=[0:0.1:10]'; ruido=ones(size(t))+0.2*randn(size(t)); % documente
figure()
plot(t,ruido)
% randn returns Normally distributed pseudorandom numbers given a size
% produce noise a step of 1 => +- 0.2*randomNumber around 1
% lsim(SYS,U,T) Simulate time response of dynamic systems to arbitrary inputs.
yruido=lsim(sysce,ruido,t); y=step(sysce,t); % encontrar respuesta a step con ruido y a step sin ruido
 
figure; subplot(211); plot(t,y); title('Respuesta sin ruido');grid; % graficar respuesta a step sin ruido

subplot(212);  plot(t,[ruido yruido]); title('Respuesta con ruido');grid % graficar respuesta a step con ruido

% se ve que el ruido no afecto casi el sistema retroalimentado



sisotool('bode',sysab)


s = tf('s');
sys1 = 10/s/(0.5*s+1);
sys2 = sys1/(s+1);
sys3 = sys1*(s+1);

