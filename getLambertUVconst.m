
function [tA_cos, A, c2, c3, si_min, si_max, si] = getLambertUVconst(ephemerisData_launch, ephemerisData_arrival, dm_plus_id,dm_minus_id)

tA_cos = (dot(ephemerisData_launch.r_vec,ephemerisData_arrival.r_vec)/(ephemerisData_launch.r_mag*ephemerisData_arrival.r_mag));

if dm_plus_id ~= 0
    A_dm_plus = 1 * sqrt(ephemerisData_launch.r_mag(dm_plus_id)*ephemerisData_arrival.r_mag(dm_plus_id)*(1+tA_cos(dm_plus_id)));
else
    A_dm_plus = 0;
end

if dm_minus_id ~= 0
    A_dm_minus = -1 * sqrt(ephemerisData_launch.r_mag(dm_minus_id)*ephemerisData_arrival.r_mag(dm_minus_id)*(1+tA_cos(dm_minus_id)));
else
    A_dm_minus = 0;
end

A_temp = [A_dm_plus,A_dm_minus];

A_id = find(A_temp~=0);
A = A_temp(A_id);

c2 = 0.5;
c3 = 1/6;
si = 0;
si_min = -4*pi;
si_max = 4*pi^2;

end