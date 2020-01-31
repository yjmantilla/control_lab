function [tm,tau] = model2points(time,signal,stableVal,p1,p2,a,b,c,d)
    [~,t1] = findClosest(signal,stableVal*p1);
    t1 = time(t1);
    [~,t2] = findClosest(signal,stableVal*p2);
    t2 = time(t2);
    
    tau = abs(a*t1 + b*t2);
    tm = abs(c*t1 + d*t2);
    
end