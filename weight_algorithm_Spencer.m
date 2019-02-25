%% Step 1 Guess W_0_old given, R, E, W_crew, W_payload
R = 1311;                                               % Range (miles)
E = 2.75;                                               % Endurance (hour)
W_crew = 340;                                           % Weight of the crew (pounds)
W_payload = 0;                                          % Weight of the payload (pounds)
W_0_old = 12100                                         % Old weight in (pounds)

%% Step 2 Find W_f/W_0
V = 943.8;                                              % Velocity (mph)
LD_max = 9;                                             % Lift max from table 3.1
LD_c = 0.866*LD_max;                                    % Lift cruise
LD_l = LD_max;                                          % Lift loiter
C_c = 0.8;                                              % Thrust Specific Fuel Consumption cruise From table 3.3
C_l = 0.7;                                              % Thrust Specific Fuel Consumption loiter From table 3.3
takeoff = 0.97;                                         % Wi/Wi-1
climb = 0.985;                                          % Wi/Wi-1
land = 0.995;                                           % Wi/Wi-1
cruise = exp(-((R*C_c)/(V*(LD_c))));                    % Wi/Wi-1
loiter = exp(-((E*C_l)/(LD_l)));                        % Wi/Wi-1
WfW0 = 1.05*(1-(takeoff*climb*land*loiter*cruise))      % 

%% Step 3 find W_e/W_0
A = 1.59;
C = -0.10;
Kvs = 1.00;
WeW0 = A * W_0_old^C * Kvs

%% Step 4

W_0_new = (W_crew + W_payload)/(1 - WfW0 - WeW0)

%% Step 5

Check = abs(W_0_new - W_0_old)/(W_0_old)