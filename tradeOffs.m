

%clear all;
%clc;

numCrew = 2;
weightPerson = 170; %pounds

weightCrew = numCrew*weightPerson;

weightPayload = 0;%4600


edurance = 2.5;
mach = 746;

W_0_range = 10000:1000:30000; %pounds
LD_range = 4:.5:14
V_c_range = linspace(0,1.7,21);
range_range = linspace(500, 2500, 21)
for i = 1:size(range_range, 2)
range = range_range(i);%Miles    
W_0 = 16028.31;

LD_max = 10;%LD_range(i);                                             % Lift max from table 3.1
LD_c = 0.866*LD_max;                                    % Lift cruise
LD_sprint = .6*LD_max;
LD_l = LD_max; % Lift loiter

V_c = V_c_range(i)*mach;
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
WeW0 = A * W_0^C * Kvs

%% Find new weight
W_0_new = (weightCrew + weightPayload)/(1 - WfW0 - WeW0)
W_0_new_range(i) = W_0_new;

change(i) = abs(W_0_new - W_0)/(W_0)

end

plot(W_0_new_range, LD_range)
