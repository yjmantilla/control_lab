diary('MantillaP06')
tic
clock

close all
s=tf('s');
w=logspace(-2,1,500)';
k = 10;
g =  k/s/(s+2.5)/(s+5);
sysab = g;
[Gm,Pm,Wcp,Wcg] = margin(g);
Gm = 20*log10(Gm);
[mag,phase]=bode(sysab,w);
figure()
bode(sysab)
mag = 20*log10(squeeze(mag));
phase = squeeze(phase);
dist    = abs(mag - (-3));
minDist = min(dist);
idx     = find(dist == minDist);
waba = w(idx);
sysce = feedback(sysab,1);
[mag2,phase2]=bode(sysce,w);
mag2 = 20*log10(squeeze(mag2));
phase2 = squeeze(phase2);
mr = max(mag2);
idr = find(mag2 == mr);
wr = w(idr);
dist    = abs(mag2 - (-3));
minDist = min(dist);
idx     = find(dist == minDist);
wabc = w(idx);
figure()
step(sysce)

close all
s=tf('s');
w=logspace(-2,2,500)';
wn = 1;
e = 0.7;
g =  wn*wn/s/(s+2*e*wn);
sysab = g;
[Gm,Pm,Wcp,Wcg] = margin(g);
Gm = 20*log10(Gm);
[mag,phase]=bode(sysab,w);
figure()
bode(sysab,w)
mag = 20*log10(squeeze(mag));
phase = squeeze(phase);
dist    = abs(mag - (-3));
minDist = min(dist);
idx     = find(dist == minDist);
waba = w(idx);
sysce = feedback(sysab,1);
[mag2,phase2]=bode(sysce,w);
figure()
bode(sysce,w)
mag2 = 20*log10(squeeze(mag2));
phase2 = squeeze(phase2);
mr = max(mag2);
idr = find(mag2 == mr);
wr = w(idr);
dist    = abs(mag2 - (-3));
minDist = min(dist);
idx     = find(dist == minDist);
wabc = w(idx);
figure()
step(sysce)





clock
toc
diary