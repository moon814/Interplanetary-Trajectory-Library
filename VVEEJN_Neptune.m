% %% constants
mu_sun = 1.32712440018e11;
mu_venus = 3.24858599e5;
mu_earth = 3.986004415e5;
mu_jupiter = 1.266865361e8;

rad_venus = 6051.8;
rad_earth = 6378.1363;
rad_jupiter = 71492;

earth_Period = 365.242189;
venus_Period = 224.6906;

% %% Set up dates
% 
% % Launch date
% 
% Y_earthBegin = 2025;
% M_earthBegin = 1; 
% D_earthBegin = 1; 
% UT_earthBegin = 0;
% Y_earthEnd = 2025; 
% M_earthEnd = 8;
% D_earthEnd = 31;
% UT_earthEnd = 0;
% 
% [earthlaunch_JD_begin] = convertDatestoJD(Y_earthBegin, M_earthBegin, D_earthBegin, UT_earthBegin);
% [earthlaunch_JD_end] = convertDatestoJD(Y_earthEnd, M_earthEnd, D_earthEnd, UT_earthEnd);
% earthLaunchDates = [earthlaunch_JD_begin:earthlaunch_JD_end];

%% PCP setup - earth to venus 
% mu_sun = 1.32712440018e11;
% 
% JD_launch_begin = earthlaunch_JD_begin; % problem 2: 2461295; 
% JD_launch_end = earthlaunch_JD_end; % problem 2: 2461415;
% JD_launch = [earthlaunch_JD_begin:earthlaunch_JD_end];
% date_launch = datetime(JD_launch,'convertfrom','juliandate');
% 
% JD_arrival_begin =  earthlaunch_JD_begin+150; %2461530; 
% JD_arrival_end = earthlaunch_JD_end+150; %2461710;
% JD_arrive = [JD_arrival_begin:JD_arrival_end];
% date_arrival = datetime(JD_arrive,'convertfrom','juliandate');
% 
% launch_T = (JD_launch - 2451545.0)/(36525);
% arrive_T = (JD_arrive - 2451545.0)/(36525);
% 
% v_mag_1 = NaN(length(launch_T),length(JD_arrive));
% C3 = NaN(length(launch_T),length(JD_arrive));
% T_0 = NaN(length(launch_T),length(JD_arrive));

%% Lambert solver
% for i = 1:length(launch_T)
%     for j = 1:length(arrive_T)
%         
%         [a_launch, e_launch, inc_launch, RAAN_launch, tA_launch, AOP_launch, M_anom_launch] = meeusEphemeris(launch_T(i), "Earth");
%         [a_arrival, e_arrival, inc_arrival, RAAN_arrival, tA_arrival, AOP_arrival, M_anom_arrival] = meeusEphemeris(arrive_T(j), "Venus");
% 
%         [r_mag_launch, r_vec_launch, v_mag_launch, v_vec_launch, dm_launch] = getStatesfromOrbParams(a_launch, e_launch, inc_launch, RAAN_launch, tA_launch, AOP_launch, M_anom_launch);
%         [r_mag_arrival, r_vec_arrival, v_mag_arrival, v_vec_arrival, dm_arrival] = getStatesfromOrbParams(a_arrival, e_arrival, inc_arrival, RAAN_arrival, tA_arrival, AOP_arrival, M_anom_arrival);
% 
%         [delta_nu(i,j), dm_plus_id, dm_minus_id] = getDeltaNuAndDM(tA_launch,tA_arrival,dm_launch,dm_arrival);
% 
%         [tA_cos, A, c2, c3, si_min, si_max, si] = getLambertUVconst_old(r_mag_launch,r_mag_arrival,r_vec_launch,r_vec_arrival,dm_plus_id,dm_minus_id);
%         
%         TOF_days = [date_arrival(j),date_launch(i)];
%         TOF_days_diff = caldiff(TOF_days,'days');
%         T_0(i,j) = abs((split(TOF_days_diff,{'days'}))*86400);
%         if T_0(i,j)<=36*86400
%             continue
%         end
% 
%         c2 = 0.5;
%         c3 = 1/6;
%         si = 0;
%         si_min = -4*pi;
%         si_max = 4*pi^2;
%         T_ls = 9e10;
%         y = 0;
%         while abs(T_ls-T_0(i,j))>1e-4
% 
%             y = r_mag_launch + r_mag_arrival + ((A*(si*c3 - 1))/sqrt(c2));
%             if y < 0 && A>0
%                 si_new = si + 0.1;
%                 y = r_mag_launch + r_mag_arrival + ((A*(si_new*c3 - 1))/sqrt(c2));
%                 c2 = (1 - cos(sqrt(si_new)))/si_new;
%                 c3 = (sqrt(si_new)-sin(sqrt(si_new)))/sqrt(si_new^3);
%             end
%             khi = sqrt(y/c2);
%             T_ls = (khi^3*c3 + A*sqrt(y))/(sqrt(mu_sun));
% 
%             if T_ls<=T_0(i,j)                       
%                 si_min = si;
%             else
%                 si_max = si;
%             end
% 
%             si = (si_max + si_min)/2;
% 
%             if (si > 1e-6)
%                 c2 = (1 - cos(sqrt(si)))/si;
%                 c3 = (sqrt(si)-sin(sqrt(si)))/sqrt(si^3);
%             elseif (si < -1e-6)
%                 c2 = (1 - cosh(sqrt(-si)))/si;
%                 c3 = (sinh(sqrt(-si))-sqrt(-si))/sqrt(-si^3);
%             else
%                 c2 = 0.5; 
%                 c3 = 1/6;
%             end
%         end
%         
%         f = 1 - (y/r_mag_launch);
%         g_dot = 1 - (y/r_mag_arrival);
%         g = A*sqrt(y/mu_sun);
% 
%         v_departure = (r_vec_arrival - (f*r_vec_launch))/g;
%         v_arrival = ((g_dot*r_vec_arrival)-r_vec_launch)/g;
% 
%         v_inf_out = v_departure - v_vec_launch;
%         v_inf_in = v_arrival - v_vec_arrival;
% 
%         v_mag_1(i,j) = norm(v_inf_in);
%         v_mag_2 = norm(v_inf_out);
%         C3(i,j) = v_mag_2^2;
%     end
% end

%% plot PCP

% date_launch_days = convertTo(date_launch,'yyyymmdd');
% date_arrival_days = convertTo(date_arrival,'yyyymmdd');
% T_0_d = T_0/86400;
% 
% figure
% hold on
% C3_contours = [2, 4, 5, 6,7,8,10, 15, 20, 30, 40, 50, 60];
% [cs1, h1]=contour(1:length(date_launch_days), 1:length(date_arrival_days), C3', C3_contours, 'r');
% clabel(cs1 ,h1);
% V_inf_contours = [3, 3.5, 4, 4.5, 5, 6, 7, 7.5, 8, 8.5,9, 11, 12, 13, 14, 15, 17,18,23, 33, 40,50, 70, 80,90];
% [cs2, h2] = contour(1:length(date_launch_days), 1:length(date_arrival_days), v_mag_1', V_inf_contours, 'b');
% clabel(cs2, h2);
% dt_contours = [10,20,30,40, 50, 60, 70, 80, 90, 100, 150,200, 250, 275,300];
% [cs3,h3] = contour(1:length(date_launch_days), 1:length(date_arrival_days), T_0_d', dt_contours, 'k');
% clabel(cs3,h3);
% title('Earth to Venus 2029-2030, Type I ');
% xlabel('Earth departure dates - 2029'); 
% ylabel('Venus arrival dates - 2030'); 
% legend('C3 (km^{2}/s^{2})','V_\infty @ Venus (km/s)','Time of Flight(days)');


%% searching algorithm 

% venus_arrival_range = [JD_arrival_begin+60:JD_arrival_begin+134]; % days when max c3 is 20 (km/s)^2
% earthLaunchDates_SA = [earthlaunch_JD_begin+60:earthlaunch_JD_begin+150];
% ega1ArrivalDates = [JD_arrival_begin+150+80:JD_arrival_begin+150+350];
% ega2launchDates = [ega1ArrivalDates(1)+(earth_Period*2):ega1ArrivalDates(end)+(earth_Period*2)];
% % ega2launchDates = [JD_arrival_begin+150+80+365+200:JD_arrival_begin+150+80+365+290];
% JOIarrivalDates = [ega2launchDates(1)+300+660:5:ega2launchDates(1)+300+1475];
earthLaunchDates_SA = [2460522.5:2460826.5];
vga1arrivalDates = [2460674.5:2460917.5]; % days when max c3 is 20 (km/s)^2
vga2launchDates = [vga1arrivalDates(1)+(venus_Period*2):vga1arrivalDates(end)+(venus_Period*2)];
ega1ArrivalDates = [vga2launchDates(1)+90:vga2launchDates(end)+600];
ega2launchDates = [ega1ArrivalDates(1)+(earth_Period*5):ega1ArrivalDates(end)+(earth_Period*5)];
% ega2launchDates = [JD_arrival_begin+150+80+365+200:JD_arrival_begin+150+80+365+290];
JupiterGADates = [ega2launchDates(1)+(earth_Period*.5):10:ega2launchDates(end)+(earth_Period*4)];
NeptuneArrivalDates = [JupiterGADates(1)+(earth_Period*4):10:JupiterGADates(end)+(earth_Period*10)];


launchPlanet = "Earth";
arrivalPlanet_VGA1 = "Venus";
launchPlanet_VGA2toEGA1 = "Venus";
arrivalPlanet_VGA2toEGA1 = "Earth";
launchPlanet_EGA2toJup = "Earth";
arrivalPlanet_EGA2toJup = "Jupiter";
launchPlanet_JuptoNep = "Jupiter";
arrivalPlanet_JuptoNep = "Neptune";

minRp_venus = (500+rad_venus);
delV_tol_VGA = 0.0005;
delV_tol_EGA = 0.09;
delV_tol = 0.45;
minRp_earth = (300+rad_earth);
minRp_jup = 200000;
max_C3 = 40;
ansCounter = 0;
valid_traj = [];
c3_id = 0;

for i = 1:length(earthLaunchDates_SA)
    
    launchDate = earthLaunchDates_SA(i);
    
    for j = 1:length(vga1arrivalDates)
        
        vga1Date = vga1arrivalDates(j);
        if (vga1Date-launchDate) < 0 || (vga1Date-launchDate) <= 90
            continue 
        end
        [delta_T_0_sec_EtoVGA1,ephemerisData_launch_EtoVGA1, ephemerisData_arrival_EtoVGA1, A_EtoVGA1] = ephemeris(launchDate, vga1Date, launchPlanet, arrivalPlanet_VGA1);
%         if delta_T_0_sec_EtoVGA1 <= 90*86400
%             continue
%         end
        [y_EtoVGA1] = lambertSolver(A_EtoVGA1,delta_T_0_sec_EtoVGA1,ephemerisData_launch_EtoVGA1.r_mag,ephemerisData_arrival_EtoVGA1.r_mag);
        [v_inf_out_earth,v_inf_in_VGA1,v_sc_departure_earth,v_sc_arrival_VGA1,C3_earth] = getVinfinityAndSome(y_EtoVGA1,ephemerisData_launch_EtoVGA1.r_mag,ephemerisData_arrival_EtoVGA1.r_mag,ephemerisData_arrival_EtoVGA1.r_vec,ephemerisData_launch_EtoVGA1.r_vec,ephemerisData_launch_EtoVGA1.v_vec,ephemerisData_arrival_EtoVGA1.v_vec,A_EtoVGA1);
        
        if C3_earth < max_C3
            
            c3_id = c3_id + 1;
            c3_dates(c3_id,1:2) = [launchDate,vga1Date];
            
%         end
%     end
% end
            vga2Date = vga2launchDates(j);
            
            for q = 1:length(ega1ArrivalDates)
                ega1ArrivalDate = ega1ArrivalDates(q);
                if (ega1ArrivalDate-vga2Date) <= 90
                    continue
                end
                [delta_T_0_sec_VGA2toEGA1,ephemerisData_launch_VGA2toEGA1,ephemerisData_arrival_VGA2toEGA1, A_VGA2toEGA1] = ephemeris(vga2Date, ega1ArrivalDate, launchPlanet_VGA2toEGA1, arrivalPlanet_VGA2toEGA1);
%                 if delta_T_0_sec_VGA2toEGA1 <= 90*86400
%                     continue
%                 end
%                 if delta_T_0_sec_VGA2toEGA1 >= 321*86400
%                     continue
%                 end
                [y_VGA2toEGA1] = lambertSolver(A_VGA2toEGA1,delta_T_0_sec_VGA2toEGA1,ephemerisData_launch_VGA2toEGA1.r_mag,ephemerisData_arrival_VGA2toEGA1.r_mag);
                [v_inf_out_VGA2,v_inf_in_EGA1,v_sc_departure_venus,v_sc_arrival_earth,C3_venus] = getVinfinityAndSome(y_VGA2toEGA1,ephemerisData_launch_VGA2toEGA1.r_mag,ephemerisData_arrival_VGA2toEGA1.r_mag,ephemerisData_arrival_VGA2toEGA1.r_vec,ephemerisData_launch_VGA2toEGA1.r_vec,ephemerisData_launch_VGA2toEGA1.v_vec,ephemerisData_arrival_VGA2toEGA1.v_vec,A_VGA2toEGA1);
                
                del_v_inf_VG = abs(norm(v_inf_in_VGA1)-norm(v_inf_out_VGA2));
                
                if del_v_inf_VG <=  delV_tol_VGA% final tolerance 0.05

                    X_VG = 2;
                    Y_VG = 1;
                    P_sc_VG = venus_Period*(X_VG/Y_VG)*86400;
                    % s/c period -> s/c SMA
                    a_sc_VG = nthroot((((P_sc_VG/(2*pi))^2)*mu_sun),3);
                    % s/c SMA, r_planet_VGA_1, avg of v_inf at VGA 1 and 2 -> v_sc_VGA_1
                    v_sc_sun_VG = sqrt(((2*mu_sun)/ephemerisData_arrival_EtoVGA1.r_mag)-(mu_sun/a_sc_VG));

                    % Solve for theta (angle between v_planet_eg1 and v_s_c)

                    % acos of v_sc_EGA_1 and avg of v_inf at EGA 1 and 2 terms
                    avg_v_inf_VG = abs(norm(v_inf_in_VGA1)+norm(v_inf_out_VGA2))/2;
                    theta_VG = acos((avg_v_inf_VG^2 + ephemerisData_arrival_EtoVGA1.v_mag^2 - v_sc_sun_VG^2)/(2*avg_v_inf_VG*ephemerisData_arrival_EtoVGA1.v_mag));

                    % calculate v_infinity_EGA_1_out and v_infinity_EGA_2_in for 0-2pi phi values 

                    % v_infinity_EGA_1_outgoing = f(theta, phi, v_infinity_EGA_1_incoming)
                    v_uV_VG = ephemerisData_arrival_EtoVGA1.v_vec/ephemerisData_arrival_EtoVGA1.v_mag;
                    N_uV_VG = (cross(ephemerisData_arrival_EtoVGA1.r_vec,ephemerisData_arrival_EtoVGA1.v_vec))/(norm(cross(ephemerisData_arrival_EtoVGA1.r_vec,ephemerisData_arrival_EtoVGA1.v_vec)));
                    c_uV_VG = cross(v_uV_VG,N_uV_VG);
                    dcm_vnc_ec_VG = [v_uV_VG, N_uV_VG, c_uV_VG];
                    phi_VG = linspace(0,2*pi,90);
                    v_inf_out_VGA1_vnc= [NaN;NaN;NaN];
                    % v_inf_out_EG1_ec = [NaN;NaN;NaN];

                    for m = 1:length(phi_VG)

                        v_inf_out_VGA1_vnc = avg_v_inf_VG*[cos(pi-theta_VG);sin(pi-theta_VG)*cos(phi_VG(m));-sin(pi-theta_VG)*sin(phi_VG(m))];

                        v_inf_out_VGA1_ec(m,:) = (dcm_vnc_ec_VG*v_inf_out_VGA1_vnc);
                        v_inf_in_VGA2(m,1:3) = v_inf_out_VGA1_ec(m,:)' + ephemerisData_arrival_EtoVGA1.v_vec - ephemerisData_launch_VGA2toEGA1.v_vec;

                        psi_vga_1 = acos((dot(v_inf_in_VGA1,v_inf_out_VGA1_ec(m,:))/(norm(v_inf_in_VGA1)*norm(v_inf_out_VGA1_ec(m,:)))));
                        psi_vga_2 = acos((dot(v_inf_out_VGA2,v_inf_in_VGA2(m,:))/(norm(v_inf_out_VGA2)*norm(v_inf_in_VGA2(m,:)))));

                        r_p_vga_1(m) = (mu_venus/norm(v_inf_out_VGA1_ec(m,:))^2)*((1/cos((pi-psi_vga_1)/2)) - 1);
                        r_p_vga_2(m) = (mu_venus/norm(v_inf_in_VGA2(m,:))^2)*((1/cos((pi-psi_vga_2)/2)) - 1);

                        r_p_vga_3(m) = (mu_venus/norm(v_inf_in_VGA2(m,:))^2)*((1/cos((pi-psi_vga_2)/2)) - 1);

                    end

                    minRp_VG1_id = find(r_p_vga_1 > minRp_venus);
                    minRp_VG2_id = find(r_p_vga_2 > minRp_venus);

                    if isempty(minRp_VG1_id)
                        continue
                    end
                    if isempty(minRp_VG2_id)
                        continue
                    end

                    minRp_vg1 = r_p_vga_1(minRp_VG1_id);
                    minRp_vg2 = r_p_vga_2(minRp_VG2_id);
                    phi_vg1 = rad2deg(phi_VG(minRp_VG1_id));
                    phi_vg2 = rad2deg(phi_VG(minRp_VG2_id));

                    phi_vg1_min = phi_vg1(1);
                    phi_vg1_max = phi_vg1(end);
                    phi_vg2_min = phi_vg2(1);
                    phi_vg2_max = phi_vg2(end);

                    % checking for phi lower bounds
                    if phi_vg1_min < phi_vg2_min
                        phi_overlap_min_VG = phi_vg2_min;
                    else
                        phi_overlap_min_VG  = phi_vg1_min;
                    end
                    if (phi_overlap_min_VG  > phi_vg1_max) || (phi_overlap_min_VG  > phi_vg2_max)
                        continue
                    end

                    % checking for phi upper bounds
                    if phi_vg1_max < phi_vg2_max
                        phi_overlap_max_VG  = phi_vg1_max;
                    else
                        phi_overlap_max_VG  = phi_vg2_max;
                    end
                    if (phi_overlap_max_VG  < phi_vg1_min) || (phi_overlap_max_VG  < phi_vg2_min)
                        continue
                    end

                    min_phi_vg1_id = find(phi_vg1 == phi_overlap_min_VG );
                    min_phi_vg2_id = find(phi_vg2 == phi_overlap_min_VG );

                    max_phi_vg1_id = find(phi_vg1 == phi_overlap_max_VG );
                    max_phi_vg2_id = find(phi_vg2 == phi_overlap_max_VG );

                    rp_vg1_range = [min_phi_vg1_id:max_phi_vg1_id];
                    rp_vg2_range = [min_phi_vg2_id:max_phi_vg2_id];

                    rp_vg1 = minRp_vg1(rp_vg1_range);
                    rp_vg2 = minRp_vg2(rp_vg2_range);

                    maxRp_vg1_id = find(minRp_vg1 == max(rp_vg1(:)));
                    maxRp_vg2_id = find(minRp_vg2 == max(rp_vg2(:)));

                    if minRp_vg1(maxRp_vg1_id) < minRp_vg2(maxRp_vg2_id)
                        optimal_phi_d_VG = phi_vg1(maxRp_vg1_id);
                        optimal_id_phi_vinf_VG = find(abs(phi_VG(:) - deg2rad(phi_vg1(maxRp_vg1_id))) < 0.0001);
                        optimal_id_rp_VG = maxRp_vg1_id;
                    else
                        optimal_phi_d_VG = phi_vg2(maxRp_vg2_id);
                        optimal_id_phi_vinf_VG = find(abs(phi_VG(:) - deg2rad(phi_vg2(maxRp_vg2_id))) < 0.0001);
                        optimal_id_rp_VG = maxRp_vg2_id;
                    end

                    opt_rp_vga_1 = r_p_vga_1(optimal_id_phi_vinf_VG);
                    opt_rp_vga_2 = r_p_vga_2(optimal_id_phi_vinf_VG);

                    if opt_rp_vga_1 >= minRp_venus                                   
                        if opt_rp_vga_2 >=  minRp_venus 
                                              
                            ega2launchDate = ega2launchDates(q);

                            for p = 1:length(JupiterGADates)

                                JupiterGADate = JupiterGADates(p);
                                [delta_T_0_sec_EGA2toJup,ephemerisData_launch_EGA2toJup,ephemerisData_arrival_EGA2toJup, A_EGA2toJup] = ephemeris(ega2launchDate, JupiterGADate, launchPlanet_EGA2toJup, arrivalPlanet_EGA2toJup);
                                if delta_T_0_sec_EGA2toJup <= 250*86400
                                    continue
                                end
                                [y_EGA2toJup] = lambertSolver(A_EGA2toJup,delta_T_0_sec_EGA2toJup,ephemerisData_launch_EGA2toJup.r_mag,ephemerisData_arrival_EGA2toJup.r_mag);
                                [v_inf_out_EGA2,v_inf_in_JGA,v_sc_departure_EGA2,v_sc_arrival_jup,C3_EGA2toJup] = getVinfinityAndSome(y_EGA2toJup,ephemerisData_launch_EGA2toJup.r_mag,ephemerisData_arrival_EGA2toJup.r_mag,ephemerisData_arrival_EGA2toJup.r_vec,ephemerisData_launch_EGA2toJup.r_vec,ephemerisData_launch_EGA2toJup.v_vec,ephemerisData_arrival_EGA2toJup.v_vec,A_EGA2toJup);

                                del_v_inf_EG_diff = abs(norm(v_inf_out_EGA2)-norm(v_inf_in_EGA1));
                                if del_v_inf_EG_diff < delV_tol_EGA
                                    
                                    for g = 1:length(NeptuneArrivalDates)
                                        NeptuneArrivalDate = NeptuneArrivalDates(g);
                                        [delta_T_0_sec_JuptoNep,ephemerisData_launch_JuptoNep,ephemerisData_arrival_JuptoNep, A_JuptoNep] = ephemeris(JupiterGADate, NeptuneArrivalDate, launchPlanet_JuptoNep, arrivalPlanet_JuptoNep);
                                        if delta_T_0_sec_JuptoNep <= 250*86400
                                            continue
                                        end
                                        [y_JuptoNep] = lambertSolver(A_JuptoNep,delta_T_0_sec_JuptoNep,ephemerisData_launch_JuptoNep.r_mag,ephemerisData_arrival_JuptoNep.r_mag);
                                        [v_inf_out_JGA,v_inf_in_NOI,v_sc_departure_EGA2,v_sc_arrival_jup,C3_JuptoNep] = getVinfinityAndSome(y_JuptoNep,ephemerisData_launch_JuptoNep.r_mag,ephemerisData_arrival_JuptoNep.r_mag,ephemerisData_arrival_JuptoNep.r_vec,ephemerisData_launch_JuptoNep.r_vec,ephemerisData_launch_JuptoNep.v_vec,ephemerisData_arrival_JuptoNep.v_vec,A_JuptoNep);

                                        if norm(v_inf_in_NOI) <= 12
                                            
                                            del_v_inf_JGA_diff = abs(norm(v_inf_out_JGA)-norm(v_inf_in_JGA));
                                            
                                            if del_v_inf_JGA_diff < delV_tol
                                            
                                                V_in_JGA = v_inf_in_JGA; %(km/s)
                                                V_in_mag_JGA = norm(V_in_JGA);
                                                V_out_JGA = v_inf_out_JGA; %(km/s)
                                                V_out_mag_JGA = norm(V_out_JGA);
                                                v_avg_JGA = (V_in_mag_JGA+V_out_mag_JGA)/2;
                                                si_JGA = acos((dot(V_in_JGA,V_out_JGA))/(V_in_mag_JGA*V_out_mag_JGA));
                                                rP_jupiter = (mu_jupiter/v_avg_JGA^2)*((1/(cos((pi-si_JGA)/2)))-1);
                                            
                                                if rP_jupiter >= minRp_jup                                                                                      
                                                    X_EG = 5;
                                                    Y_EG = 1;
                                                    P_sc_EG = earth_Period*(X_EG/Y_EG)*86400;
                                                    % s/c period -> s/c SMA
                                                    a_sc_EG = nthroot((((P_sc_EG/(2*pi))^2)*mu_sun),3);
                                                    % s/c SMA, r_planet_EGA_1, avg of v_inf at EGA 1 and 2 -> v_sc_EGA_1
                                                    v_sc_sun_EG = sqrt(((2*mu_sun)/ephemerisData_arrival_VGA2toEGA1.r_mag)-(mu_sun/a_sc_EG));

                                                    % Solve for theta (angle between v_planet_eg1 and v_s_c)

                                                    % acos of v_sc_EGA_1 and avg of v_inf at EGA 1 and 2 terms
                                                    avg_v_inf_EG = (norm(v_inf_in_EGA1) + norm(v_inf_out_EGA2))/2;
                                                    theta_EG = acos((avg_v_inf_EG^2 + ephemerisData_arrival_VGA2toEGA1.v_mag^2 - v_sc_sun_EG^2)/(2*avg_v_inf_EG*ephemerisData_arrival_VGA2toEGA1.v_mag));

                                                    % calculate v_infinity_EGA_1_out and v_infinity_EGA_2_in for 0-2pi phi values 

                                                    % v_infinity_EGA_1_outgoing = f(theta, phi, v_infinity_EGA_1_incoming)
                                                    v_uV_EG = ephemerisData_arrival_VGA2toEGA1.v_vec/ephemerisData_arrival_VGA2toEGA1.v_mag;
                                                    N_uV_EG = (cross(ephemerisData_arrival_VGA2toEGA1.r_vec,ephemerisData_arrival_VGA2toEGA1.v_vec))/(norm(cross(ephemerisData_arrival_VGA2toEGA1.r_vec,ephemerisData_arrival_VGA2toEGA1.v_vec)));
                                                    c_uV_EG = cross(v_uV_EG,N_uV_EG);
                                                    dcm_vnc_ec_EG = [v_uV_EG, N_uV_EG, c_uV_EG];
                                                    phi_EG = linspace(0,2*pi,90);
                                                    v_inf_out_EGA1_vnc= [NaN;NaN;NaN];
                                                    % v_inf_out_EG1_ec = [NaN;NaN;NaN];

                                                    for w = 1:length(phi_EG)

                                                        v_inf_out_EGA1_vnc = avg_v_inf_EG*[cos(pi-theta_EG);sin(pi-theta_EG)*cos(phi_EG(w));-sin(pi-theta_EG)*sin(phi_EG(w))];

                                                        v_inf_out_EGA1_ec(w,:) = (dcm_vnc_ec_EG*v_inf_out_EGA1_vnc);
                                                        v_inf_in_EGA2(w,1:3) = v_inf_out_EGA1_ec(w,:)' + ephemerisData_arrival_VGA2toEGA1.v_vec - ephemerisData_launch_EGA2toJup.v_vec;

                                                        psi_ega_1 = acos((dot(v_inf_in_EGA1,v_inf_out_EGA1_ec(w,:))/(norm(v_inf_in_EGA1)*norm(v_inf_out_EGA1_ec(w,:)))));
                                                        psi_ega_2 = acos((dot(v_inf_out_EGA2,v_inf_in_EGA2(w,:))/(norm(v_inf_out_EGA2)*norm(v_inf_in_EGA2(w,:)))));

                                                        r_p_ega_1(w) = (mu_earth/norm(v_inf_out_EGA1_ec(w,:))^2)*((1/cos((pi-psi_ega_1)/2)) - 1);
                                                        r_p_ega_2(w) = (mu_earth/norm(v_inf_in_EGA2(w,:))^2)*((1/cos((pi-psi_ega_2)/2)) - 1);

                                                    end

                                                    minRp_EG1_id = find(r_p_ega_1 > minRp_earth);
                                                    minRp_EG2_id = find(r_p_ega_2 > minRp_earth);

                                                    if isempty(minRp_EG1_id)
                                                        continue
                                                    end
                                                    if isempty(minRp_EG2_id)
                                                        continue
                                                    end

                                                    minRp_eg1 = r_p_ega_1(minRp_EG1_id);
                                                    minRp_eg2 = r_p_ega_2(minRp_EG2_id);
                                                    phi_eg1 = rad2deg(phi_EG(minRp_EG1_id));
                                                    phi_eg2 = rad2deg(phi_EG(minRp_EG2_id));

                                                    phi_eg1_min = phi_eg1(1);
                                                    phi_eg1_max = phi_eg1(end);
                                                    phi_eg2_min = phi_eg2(1);
                                                    phi_eg2_max = phi_eg2(end);

                                                    % checking for phi lower bounds
                                                    if phi_eg1_min < phi_eg2_min
                                                        phi_overlap_min_EG = phi_eg2_min;
                                                    else
                                                        phi_overlap_min_EG = phi_eg1_min;
                                                    end
                                                    if (phi_overlap_min_EG > phi_eg1_max) || (phi_overlap_min_EG > phi_eg2_max)
                                                        continue
                                                    end

                                                    % checking for phi upper bounds
                                                    if phi_eg1_max < phi_eg2_max
                                                        phi_overlap_max_EG = phi_eg1_max;
                                                    else
                                                        phi_overlap_max_EG = phi_eg2_max;
                                                    end
                                                    if (phi_overlap_max_EG < phi_eg1_min) || (phi_overlap_max_EG < phi_eg2_min)
                                                        continue
                                                    end

                                                    min_phi_eg1_id = find(phi_eg1 == phi_overlap_min_EG);
                                                    min_phi_eg2_id = find(phi_eg2 == phi_overlap_min_EG);

                                                    max_phi_eg1_id = find(phi_eg1 == phi_overlap_max_EG);
                                                    max_phi_eg2_id = find(phi_eg2 == phi_overlap_max_EG);

                                                    rp_eg1_range = [min_phi_eg1_id:max_phi_eg1_id];
                                                    rp_eg2_range = [min_phi_eg2_id:max_phi_eg2_id];

                                                    rp_eg1 = minRp_eg1(rp_eg1_range);
                                                    rp_eg2 = minRp_eg2(rp_eg2_range);

                                                    [maxRp_eg1_id] = find(minRp_eg1 == max(rp_eg1(:)));
                                                    [maxRp_eg2_id] = find(minRp_eg2 == max(rp_eg2(:)));

                                                    if minRp_eg1(maxRp_eg1_id) < minRp_eg2(maxRp_eg2_id)
                                                        optimal_phi_d_EG = phi_eg1(maxRp_eg1_id);
                                                        optimal_id_phi_vinf_EG = find(abs(phi_EG(:) - deg2rad(phi_eg1(maxRp_eg1_id))) < 0.0001);
                                                        optimal_id_rp_EG = maxRp_eg1_id;
                                                    else
                                                        optimal_phi_d_EG = phi_eg2(maxRp_eg2_id);
                                                        optimal_id_phi_vinf_EG = find(abs(phi_EG(:) - deg2rad(phi_eg2(maxRp_eg2_id))) < 0.0001);
                                                        optimal_id_rp_EG = maxRp_eg2_id;
                                                    end

                                                    v_inf_out_EGA1_final = v_inf_out_EGA1_ec(optimal_id_phi_vinf_EG,:);
                                                    v_inf_in_EGA2_final = v_inf_in_EGA2(optimal_id_phi_vinf_EG,:);

                                                    opt_rp_ega_1 = r_p_ega_1(optimal_id_phi_vinf_EG);
                                                    opt_rp_ega_2 = r_p_ega_2(optimal_id_phi_vinf_EG);

                                                    if opt_rp_ega_1 >= minRp_earth                                   
                                                        if opt_rp_ega_2 >=  minRp_earth                           
                                                            ansCounter = ansCounter+1;
                                                            valid_traj(ansCounter,1:7) = [launchDate,vga1Date,vga2Date,ega1ArrivalDate,ega2launchDate,JupiterGADate,NeptuneArrivalDate];
                                                            C3_valid(ansCounter,:) = C3_earth;
                                                            v_inf_VG_diff(ansCounter,:)  = del_v_inf_VG;  
                                                            v_inf_JG_diff(ansCounter,:) = del_v_inf_JGA_diff;
                                                            v_inf_EG_diff(ansCounter,:) = del_v_inf_EG_diff;
                                                            valid_phi_VG(ansCounter,:) = optimal_phi_d_VG;
                                                            valid_phi_EG(ansCounter,:) = optimal_phi_d_EG;
                                                            valid_v_inf_NOI(ansCounter,:) = v_inf_in_NOI;                                                          
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end   
                   end
               end
            end
        end
    end
end




valid_traj_dates = datetime(valid_traj,'convertfrom','juliandate');
% min_C3 = find(C3_valid==min(C3_valid(:)));
% min_vin_inf_EG_diff_id = find(v_inf_EG_diff==min(v_inf_EG_diff(:)));
% min_vin_inf_EG_diff = v_inf_EG_diff(min_vin_inf_EG_diff_id);
% min_vin_inf_venus_diff_id = find(v_inf_venus_diff==min(v_inf_venus_diff(:)));
% min_vin_inf_venus_diff = v_inf_venus_diff(min_vin_inf_venus_diff_id);
% min_C3_id = find(C3_valid==min(C3_valid(:)));
% min_C3 = C3_valid(min_C3_id);
for k = 1:length(valid_v_inf_NOI)
    v_inf_NOI_mag(k,:) = norm(valid_v_inf_NOI(k,:));
end