
function [delta_nu,dm_plus_id,dm_minus_id] = getDeltaNuAndDM(tA_launch,tA_arrival,dm_launch,dm_arrival)

delta_nu = dm_arrival-dm_launch;
% delta_nu = wrapTo2Pi(tA_arrival)-wrapTo2Pi(tA_launch);

if delta_nu < 0
    delta_nu = (2*pi)+ delta_nu;
end

dm_plus_id = find(rad2deg(delta_nu)<180);
dm_minus_id = find(rad2deg(delta_nu)>180);
    
end