% lambert solver
function [y] = lambertSolver(A,delta_T_0_sec,r_mag_launch,r_mag_arrival)

% % Lambert solver 
mu_sun = 1.32712440018e11;

T_0 = delta_T_0_sec;
c2 = 0.5;
c3 = 1/6;
si = 0;
si_min = -4*pi;
si_max = 4*pi^2;
T_ls = 9e10;
y = 0;
iter = 1;

while abs(T_ls-T_0)>1e-5

    y = r_mag_launch + r_mag_arrival + ((A*(si*c3 - 1))/sqrt(c2));           
    if y < 0 && A>0
        si_min = si_min + 0.1;
        si = (si_max + si_min)/2;        
        if (si > 1e-6)
            c2 = (1 - cos(sqrt(si)))/si;
            c3 = (sqrt(si)-sin(sqrt(si)))/sqrt(si^3);
        elseif (si < -1e-6)
            c2 = (1 - cosh(sqrt(-si)))/si;
            c3 = (sinh(sqrt(-si))-sqrt(-si))/sqrt(-si^3);
        else
            c2 = 0.5; 
            c3 = 1/6;
        end
        y = r_mag_launch + r_mag_arrival + ((A*(si*c3 - 1))/sqrt(c2));
    end
    
    khi = sqrt(y/c2);
    T_ls = (khi^3*c3 + A*sqrt(y))/(sqrt(mu_sun));

    if T_ls<=T_0                         
        si_min = si;
    else
        si_max = si;
    end

    si = (si_max + si_min)/2;

    if (si > 1e-6)
        c2 = (1 - cos(sqrt(si)))/si;
        c3 = (sqrt(si)-sin(sqrt(si)))/sqrt(si^3);
    elseif (si < -1e-6)
        c2 = (1 - cosh(sqrt(-si)))/si;
        c3 = (sinh(sqrt(-si))-sqrt(-si))/sqrt(-si^3);
    else
        c2 = 0.5; 
        c3 = 1/6;
    end
    
    iter = iter + 1; 
    if iter > 500 
        break 
    end 
end

end