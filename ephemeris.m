% ephemeris function
function [delta_T_0_sec, ephemerisData_launch, ephemerisData_arrival, A] = ephemeris(launchDate_JD, arrivalDate_JD, launchPlanet, arrivalPlanet)

% date_launch = datetime(launchDate_JD,'convertfrom','juliandate');
% date_arrival = datetime(arrivalDate_JD,'convertfrom','juliandate');
% dates = [date_launch, date_arrival];
delta_T_0_days = abs(arrivalDate_JD-launchDate_JD);
[delta_T_0_sec] = delta_T_0_days*86400;

% [delta_T_0_sec] = (split(delta_T_0_days,{'days'}))*86400;

launch_T = (launchDate_JD - 2451545.0)/(36525);
arrive_T = (arrivalDate_JD - 2451545.0)/(36525);

[a_launch, e_launch, inc_launch, RAAN_launch, tA_launch, AOP_launch, M_anom_launch] = meeusEphemeris(launch_T, launchPlanet);
[a_arrival, e_arrival, inc_arrival, RAAN_arrival, tA_arrival, AOP_arrival, M_anom_arrival] = meeusEphemeris(arrive_T, arrivalPlanet);

[ephemerisData_launch.r_mag, ephemerisData_launch.r_vec, ...
 ephemerisData_launch.v_mag, ephemerisData_launch.v_vec, dm_launch] ...
    = getStatesfromOrbParams(a_launch, e_launch, inc_launch, RAAN_launch, tA_launch, AOP_launch, M_anom_launch);
[ephemerisData_arrival.r_mag, ephemerisData_arrival.r_vec, ...
 ephemerisData_arrival.v_mag, ephemerisData_arrival.v_vec, dm_arrival] ...
    = getStatesfromOrbParams(a_arrival, e_arrival, inc_arrival, RAAN_arrival, tA_arrival, AOP_arrival, M_anom_arrival);

[delta_nu, dm_plus_id, dm_minus_id] = getDeltaNuAndDM(tA_launch,tA_arrival,dm_launch,dm_arrival);

[tA_cos, A, c2, c3, si_min, si_max, si] = getLambertUVconst(ephemerisData_launch, ephemerisData_arrival, dm_plus_id, dm_minus_id);
 
end
