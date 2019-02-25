clear variables; close all; clear windows; clc;
%%%Inputs
syms w_0_guess
converge_goal = 0.1; %%percent
n_crew = 2;
n_passenger = 0;
w_missile = 0; %%lbs
w_payload_additional = 0; %%lbs
%A = 1.65;
A = 1.59;
%C = -0.077;
C = -0.10;
K_vs = 1;
range = 1500; % 2000;%
endurance_loiter = 2.5; %%hrs
v_cruise = .9*1125; % 350;%
LD_max = 10;
%Allison T56 Series IV
% FC = 2412; %%lbs fuel/hr
% P = 5100; %%equivallent shaft hp
% FOS_fuel = 1.06;
%%%%Calculations%%%%
w_payload = w_missile + w_payload_additional;
%%empty weight
empty_weight_ratio = A * K_vs * (w_0_guess^C);
%%SFC/Flight Conditions
%C_power = FC/P;
%SFC = C_power * v_cruise * 1.4666667 / 0.9;
C_BHP = FC/P;
SFC = C_BHP * v_cruise * 1.4666667 / (550 * 0.9); %%lbs-fuel/hr / lbs-thrust
LD_cruise = LD_max;
LD_loiter = 0.866*LD_max;
%%Mission Legs
takeoff_ratio = 0.97;
climb_ratio = 0.985;
cruise_ratio = exp( -range*SFC/(v_cruise * LD_cruise) );
loiter_ratio = exp( -endurance_loiter*SFC/(LD_loiter) );
landing_ratio = 0.995;
%%Totals
fuel_ratio = FOS_fuel*(1 - takeoff_ratio*climb_ratio*cruise_ratio*loiter_ratio*landing_ratio);
w_0 = ( (n_crew + n_passenger)*170 + w_payload ) / ( 1 - fuel_ratio - empty_weight_ratio );
%%lbs
difference = (w_0 - w_0_guess)*100/w_0_guess; %% percent
figure();
fplot([w_0, w_0_guess], [90000 110000]) ;
%fplot([w_0, w_0_guess], [110000 130000]) ;
title("W_0 Calculated");
xlabel("Takeoff Weight Guessed (lbs)");
ylabel("Takeoff Weight Calculated (lbs)");
legend("W_0 Calculated", "W_0 Guessed");
% figure();
% fplot([difference], [10000 120000]) ;
% title("Difference between Estimated and Calculated W_0");
% xlabel("lbs guessed");
% ylabel("percent");
w_0_calculated = vpasolve(difference == converge_goal, w_0_guess);
answer = eval(w_0_calculated);
final_empty_weight_ratio = A * K_vs * answer^C;
fprintf("Takeoff Weight: " + answer + " lbs \n")
fprintf("Empty Weight Ratio: " + final_empty_weight_ratio + " \n")
fprintf("Empty Weight: " + final_empty_weight_ratio*answer + " lbs \n")
fprintf("Fuel Weight Ratio: " + fuel_ratio + " \n")
fprintf("Fuel Weight: " + fuel_ratio*answer + " l")