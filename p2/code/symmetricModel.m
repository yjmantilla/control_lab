function [tau1,tau2,tm,x] = symmetricModel(numberOfTests,plant_kp,plant_in,plant_out,plant_stableIndex,time,start)
% start 0 < x < 0.5
% start 1e-3 works, realmin does not, nor 1e-10
x = linspace(start,0.5-start,numberOfTests);
errors = zeros(length(x),1);
s = tf('s');
taus1 = zeros(length(x),1);
taus2 = zeros(length(x),1);
tms = zeros(length(x),1);
for i = 1:length(x)
[~,tx] = findClosest(plant_out,plant_out(plant_stableIndex)*x(i)); 
tx = time(tx);
[~,t50] = findClosest(plant_out,plant_out(plant_stableIndex)*0.5); 
t50 = time(t50);
[~,t1mx] = findClosest(plant_out,plant_out(plant_stableIndex)*(1-x(i))); 
t1mx = time(t1mx);
a = (-0.6240*tx + 0.9866*t50 - 0.3626*t1mx)/(0.3533*tx - 0.7036*t50 + 0.3503*t1mx);
a = abs(a);
tauPP = (t1mx-tx)/(0.9866+0.7036*a);
tauPP = abs(tauPP);
tau1 = tauPP;
tau2 = a*tauPP;
tm = t1mx - (1.3421+1.3455*a)*tauPP;
if tm < 0
    tm = 0;
end
sys = (plant_kp * exp(-tm*s))/((tau1*s+1)*(tau2*s+1));
sys_out = step(sys,time)';
errors(i) = immse(sys_out,plant_out/max(plant_in));
taus1(i) = tau1;
taus2(i) = tau2;
tms(i) = tm;
end

[~,index] = min(errors);

tau1 = taus1(index);
tau2 = taus2(index);
tm = tms(index);
x = x(index);

end