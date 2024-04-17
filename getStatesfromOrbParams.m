
function [r_mag, r_vec, v_mag, v_vec, dm] = getStatesfromOrbParams(a, e, inc, RAAN, tA, AOP, M_anom)

theta = AOP + tA; 

C_1 = cos(RAAN)*cos(theta) - sin(RAAN)*cos(inc)*sin(theta);
C_2 = -cos(RAAN)*sin(theta) - sin(RAAN)*cos(inc)*cos(theta);
C_3 = sin(RAAN)*sin(inc);
C_4 = sin(RAAN)*cos(theta) + cos(RAAN)*cos(inc)*sin(theta);
C_5 = -sin(RAAN)*sin(theta) + cos(RAAN)*cos(inc)*cos(theta);
C_6 = -cos(RAAN)*sin(inc);
C_7 = sin(inc)*sin(theta);
C_8 = sin(inc)*cos(theta);
C_9 = cos(inc);

p = a*(1-e^2);
r_mag = (p)/(1+e*cos(tA));
% mu = 1.32712440018e11; %-- use for 6008 (earth-sun system)
mu = 3.986004415e5; % -- using for 6080 

h = sqrt(mu*p);

r_1 = r_mag*C_1;
r_2 = r_mag*C_4;
r_3 = r_mag*C_7;
r_vec = [r_1;r_2;r_3];
r_mag = norm(r_vec);

v_1 = (r_1*h*e*sin(tA)/(r_mag*p)) + (h/r_mag)*C_2;
v_2 = (r_2*h*e*sin(tA)/(r_mag*p)) + (h/r_mag)*C_5;
v_3 = (r_3*h*e*sin(tA)/(r_mag*p)) + (h/r_mag)*C_8;
v_vec = [v_1;v_2;v_3];
v_mag = norm(v_vec);

dm = atan2(r_2,r_1);
% if dm < 0
%     dm = (2*pi)-abs(dm);
% end

end