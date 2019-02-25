

clear all;
clc;
close all;

syms w_0_guess
converge_goal = 0.1; %%percent
numCrew = 2;
weightPerson = 170; %pounds

weightCrew = numCrew*weightPerson;

weightPayload = 0;%4600

range = 1500;%Miles
edurance = 2.5;
mach = 746;

LD_max = 10;                                             % Lift max from table 3.1
LD_c = 0.866*LD_max;                                    % Lift cruise
LD_sprint = .6*LD_max;
LD_l = LD_max; % Lift loiter

V_c = .9*mach;
V_sprint = 1.5*mach;
C_c = 0.8;                                              % Thrust Specific Fuel Consumption cruise From table 3.3
C_sprint = 1.1;
C_l = 0.7;                                              % Thrust Specific Fuel Consumption loiter From table 3.3
takeoff = 0.97;                                         % Wi/Wi-1
climb = 0.985;                                          % Wi/Wi-1
land = 0.995;                                           % Wi/Wi-1
cruise = exp(-(((.75*range)*C_c)/(V_c*(LD_c))));                    % Wi/Wi-1
sprint = exp(-(((.25*range)*C_sprint)/(V_sprint*(LD_sprint)))); 
loiter = exp(-((edurance*C_l)/(LD_l)));                        % Wi/Wi-1
WfW0 = 1.05*(1-(takeoff*climb*land*loiter*cruise*sprint))      % 


%% W_e/W_0
A = 1.59;
C = -0.10;
Kvs = 1.00;
WeW0 = A * w_0_guess^C * Kvs

%% Find new weight
w_0 = (weightCrew + weightPayload)/(1 - WfW0 - WeW0)

difference = (w_0 - w_0_guess)*100/w_0_guess; %% percent
figure();
fplot([w_0, w_0_guess], [15000 20000]) ;
%fplot([w_0, w_0_guess], [110000 130000]) ;
title("W_0 Calculated");
xlabel("Takeoff Weight Guessed (lbs)");
ylabel("Takeoff Weight Calculated (lbs)");
legend("W_0 Calculated", "W_0 Guessed");

w_0_calculated = vpasolve(difference == converge_goal, w_0_guess);
answer = eval(w_0_calculated);
final_empty_weight_ratio = A * Kvs * answer^C;
fprintf("Takeoff Weight: " + answer + " lbs \n")
fprintf("Empty Weight Ratio: " + final_empty_weight_ratio + " \n")
fprintf("Empty Weight: " + final_empty_weight_ratio*answer + " lbs \n")
fprintf("Fuel Weight Ratio: " + WfW0 + " \n")
fprintf("Fuel Weight: " + WfW0*answer + " l")