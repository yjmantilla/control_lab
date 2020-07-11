function [p,pi,pid] = zn_design(ku,tu)
    p.kc = ku/2;
    p.bp = 100/p.kc;
    
    pi.kc = ku/2.2;
    pi.bp = 100/pi.kc;
    pi.ti = tu/1.2;
    
    pid.kc = ku/1.7;
    pid.bp = 100/pid.kc;
    pid.ti = tu/2;
    pid.td = tu/8;
    
end
