
function [a, e, inc, RAAN, tA, AOP, M_anom] = meeusEphemeris(T, planet)

[L_coeff, a_coeff, e_coeff,inc_coeff,RAAN_coeff,LOP_coeff] = getCoeffs(planet);


L = L_coeff(1) + L_coeff(2)*T + L_coeff(3)*T^2 + L_coeff(4)*T^3;

a = a_coeff(1) + a_coeff(2)*T + a_coeff(3)*T^2 + a_coeff(4)*T^3;

e = e_coeff(1) + e_coeff(2)*T + e_coeff(3)*T^2 + e_coeff(4)*T^3;

inc = inc_coeff(1) + inc_coeff(2)*T + inc_coeff(3)*T^2 + inc_coeff(4)*T^3;

RAAN = RAAN_coeff(1) + RAAN_coeff(2)*T + RAAN_coeff(3)*T^2 + RAAN_coeff(4)*T^3;

LOP = LOP_coeff(1) + LOP_coeff(2)*T + LOP_coeff(3)*T^2 + LOP_coeff(4)*T^3;

AOP = LOP - RAAN;

M_anom = L - LOP;

C_cen = (((2*e) - (e^3/4) + ((5*e^5)/96)))*sin(M_anom) + (((5*e^2)/4) - ((11*e^4)/24))*sin(2*M_anom) + (((13*e^3)/12) - ((43*e^5)/64))*sin(3*M_anom) + ((103*e^4*sin(4*M_anom))/96) + ((1097*e^5*sin(5*M_anom))/960);

tA = M_anom + C_cen;

end