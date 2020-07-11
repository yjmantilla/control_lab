
%plot(ScopeData.time,ScopeData.signals.values)
Tu_s = 330 -290;
Ku_s = 8.29;

Tu_p = 38;
Ku_p = 10;
[p_s,pi_s,pid_s] = zn_design(Ku_s,Tu_s);
[p_p,pi_p,pid_p] = zn_design(Ku_p,Tu_p);

step_practica_alfaro_p = stepinfo(resp_practica_alfaro_p.signals.values,resp_practica_alfaro_p.time);
step_practica_alfaro_pi = stepinfo(resp_practica_alfaro_pi.signals.values,resp_practica_alfaro_pi.time);
step_practica_alfaro_pid = stepinfo(resp_practica_alfaro_pid.signals.values,resp_practica_alfaro_pid.time);

step_sim_alfaro_p = stepinfo(resp_sim_alfaro_p.signals.values,resp_sim_alfaro_p.time);
step_sim_alfaro_pi = stepinfo(resp_sim_alfaro_pi.signals.values,resp_sim_alfaro_pi.time);
step_sim_alfaro_pid = stepinfo(resp_sim_alfaro_pid.signals.values,resp_sim_alfaro_pid.time);

step_practica_zn_p = stepinfo(resp_practica_zn_p.signals.values,resp_practica_zn_p.time);
step_practica_zn_pi = stepinfo(resp_practica_zn_pi.signals.values,resp_practica_zn_pi.time);
step_practica_zn_pid = stepinfo(resp_practica_zn_pid.signals.values,resp_practica_zn_pid.time);

step_sim_zn_p = stepinfo(resp_sim_zn_p.signals.values,resp_sim_zn_p.time);
step_sim_zn_pi = stepinfo(resp_sim_zn_pi.signals.values,resp_sim_zn_pi.time);
step_sim_zn_pid = stepinfo(resp_sim_zn_pid.signals.values,resp_sim_zn_pid.time);

step_exp_zn_p = stepinfo(resp_exp_zn_p.signals.values,resp_exp_zn_p.time);
step_exp_zn_pi = stepinfo(resp_exp_zn_pi.signals.values,resp_exp_zn_pi.time);
step_exp_zn_pid = stepinfo(resp_exp_zn_pid.signals.values,resp_exp_zn_pid.time);
save zn.mat
