close all
clear all
addpath('C:\Users\user\Desktop\code\control_lab\')
path = 'C:\Users\user\Desktop\code\control_lab\p5\';
cd(path)

wn = 1;

s = tf('s');
figure(1)
for damp= 0:0.1:1.2
sys = (wn^2)/(s*(s+2*damp*wn));
[r,k] = rlocus(sys);
hold on
rlocus(sys);
end
figure(2)
for damp= 0:0.1:1.2
sys = (wn^2)/(s*(s+2*damp*wn));
sys2 = sys/(1+sys);
hold on
pzmap(sys2);
[r,k] = rlocus(sys);
%rlocus(sys);
end


figure(3)
%for a = 0:0.2:3
a=3;
sys = s^2 + a*s +1;
[r,k] = rlocus(sys);
hold on
rlocus(sys);
%end

figure(4)

sys = (1)/(s*(s+2));
rlocus(sys);

figure(5)
sys = sys/(1+s/5);
rlocus(sys);

figure(6)
sys = sys*(1+s/3);
rlocus(sys);
