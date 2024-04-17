% function calculates the necessary departure and arrival velocities at
% specified planets state (position, velocity) and other quantities like y
% and A
function [v_inf_out,v_inf_in,v_sc_departure,v_sc_arrival,C3] = getVinfinityAndSome(y,r_mag_launch,r_mag_arrival,r_vec_arrival,r_vec_launch,v_vec_launch,v_vec_arrival,A)

mu_sun = 1.32712440018e11;

f = 1 - (y/r_mag_launch);
g_dot = 1 - (y/r_mag_arrival);
g = A*sqrt(y/mu_sun);

v_sc_departure = (r_vec_arrival - (f*r_vec_launch))/g;
v_sc_arrival = ((g_dot*r_vec_arrival)-r_vec_launch)/g;

v_inf_out = v_sc_departure - v_vec_launch;
v_inf_in = v_sc_arrival - v_vec_arrival;

v_mag_1 = norm(v_inf_in);
v_mag_2 = norm(v_inf_out);

C3 = v_mag_2^2;

% found from previous earth ot venus transfer
% v_inf_in_venus = [3.64174993634358;-3.84270636476666;-3.26879817785136];
% v_inf_out_venus = [2.412335895892081;-5.733312715017199;-0.172982521589608];
% v_inf_in_venus_mag = norm(v_inf_in_venus);
% v_inf_out_venus_mag = norm(v_inf_out_venus);
% 
% si_hyp = acos((dot(v_inf_in_venus,v_inf_out_venus))/(v_inf_in_venus_mag*v_inf_out_venus_mag));
% 
% r_p = (mu_arrival/v_inf_in_venus_mag^2)*((1/(cos((pi-si_hyp)/2)))-1)
% r_arrival = 6051.8;
% r_p_alt = r_p - r_arrival
% TOF_days = T_0/86400;

% energy_sc = ((norm(v_arrival))^2/2) - (mu_sun/(r_mag_arrival))
% temp1(i) = ((v_inf_sq*rp(i))/mu_venus);
%     ro(i) = acos(1/(1+temp1(i)));
%     psi(i) = pi - 2*ro(i);
% [12818.9705105576] km  altitude 
% energy_sc = ((norm(v_departure))^2/2) - (mu_sun/(r_mag_launch))
end
