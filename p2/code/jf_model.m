function [jf_out,jf_sys,jf_tau] = jf_model(time,plant_kp,jf_tm,t70,t90)
s = tf('s');
jf_nabla = (t90-t70)/(t90-jf_tm);
if jf_nabla <= 0.4771
    jf_damp = sqrt((0.4844651-0.75323499*jf_nabla)/(1-2.0946464*jf_nabla));
else
    jf_damp = 13.9352;
end
jf_tau = (t90-jf_tm)/(-0.424301+4.62533*jf_damp-2.65412*exp(-jf_damp));
jf_sys = (plant_kp * exp(-jf_tm*s))/(((jf_tau)^2)*(s^2)+2*jf_damp*jf_tau*s+1);
jf_out = step(jf_sys,time)';
end