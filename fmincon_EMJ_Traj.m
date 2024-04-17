clear all
clc
close all

global MU AU
MU = 1.32712440018e11;
AU = 1.49597870700e8;
day2sec  = 86400;

%Setting up the bounds for Launch, MGA, and JOI
Launch0 = 2462868;
LaunchF = Launch0 + 100;

%Mars Gravity Assist Window
MGA0 = Launch0 + 50;
MGAF = LaunchF + 200;

%Jupiter Orbit Insertion Window
JOI0 = MGA0 + 500;
JOIF = MGAF + 1000;

%Initial State Vector selected to be exactly in the middle of the allowable
%range.
X = [mean([Launch0,LaunchF]), mean([MGA0, MGAF]), mean([JOI0, JOIF])];

%This is the cost function used within fmincon. It tries to minimize a
%weighted combination of C3, Arrival-Vinfinity at Jupiter, the difference
%in V-infinity at MGA, and the TOF. 

%Cost Function = C3*w1 + diff_MGA*w2 + w3*VJOI_mag + TOF*w4;

%Weights  
%W1 is the weight for C3
%W2 is the weight for Delta_V_inf @ MGA
%W3 is the weight for V_inf_arrival @ JOI
%W4 is the weight for TOF in years. 
W = [1, 30, 1, 1];

%Nominal Evaluation
f2 = @(X)EMJ_traj(X, W);

%Options
A = [];
b = [];
Aeq = [];
beq = [];

%You can set the bounds to be the searching range, or set them to be open
%(unlimited)

% lb  = [Launch0, MGA0, JOI0];
% ub  = [LaunchF, MGAF, JOIF];

lb = [];
ub = [];

%fmincon options. See help fminopts for additional help. 
fminopts = optimoptions(@fmincon,'TolCon', 1e-12, 'Display', 'iter',...
    'FinDiffType', 'central', 'TolFun', 1e-12, 'TolX', 1e-12);

%Running fmincon.
[Xout, fval] = fmincon(f2, X, A, b, Aeq, beq, lb, ub, [], fminopts);

%Optimal Dates from fmincon
Dates = [Xout];

%What are the parameters of the optimal trajectory that fmicncon found
%based on the weights, date ranges, etc?
[C3, diff_MGA, VJOI_mag, alt_Mars, TOFy] = EMJ_traj_output(Xout);

fprintf('Trajectory Performance\n')
fprintf('**************************************\n')
fprintf('C3 = %3.3f km^2/s^2\n', C3)
fprintf('V-infinity difference at MGA = %2.3f km/s\n', diff_MGA);
fprintf('Altitude at MGA = %2.3f km/s\n', alt_Mars);
fprintf('V-infinity arrival at Jupiter = %2.3f km/s\n', VJOI_mag);
fprintf('TOF = %2.3f years\n', TOFy)

































