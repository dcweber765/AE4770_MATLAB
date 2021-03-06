clear all;

h = 30000;
w = 1.6028e+04;
S_wing = w/50;
alpha = 2;
alpha_0 = -2.35;

g = 32.15;                          %Gravitational acceleration [ft/sec^2]
R = 3089.2;                         %Specific gas constant for dry air [lb*ft/sl*K]
a1 = -1.9812*10^-3;                 %Temperature gradient [K/ft]
rho_sl = 0.0023769;                 %Rho at sea level [sl/ft^3]
tau_sl = 288.16;                    %Standard temperature at sea level [K]



if h<=36089.2                       %For altitudes 0 - 36089.2 km
    tau = tau_sl + a1*h;
    rho = rho_sl * (tau/tau_sl).^(-1-g/(a1*R));
elseif h>36089.2                    %For altitude 36089.2 - 82021 km
    tau_11 = tau_sl + a1*36089.2;   %Temperature at 36089.2 km
    rho_11 = rho_sl * (tau/tau_sl).^(-1-g/(a1*R));    %Pressure at 36089.2 km
    
    tau = tau_11;
    rho = rho_11*exp((-g*(h-36089.2))/(R*tau_11));
end

ws_TO = 50; 
Mach = .9;
v_cruise = Mach *1125;%ft/s
C_l_alpha = 2*pi;%From XFLR5 ms313
LAMDA = 40;
e = .85;

AR = 3.185;
Beta = sqrt(1-Mach^2);
eta = C_l_alpha/(2*pi/Beta);

C_L_max = .9*1.3*cosd(44.26); %add flap 

C_L_alpha = (2*pi*AR)/(2+sqrt(4+((AR^2*Beta^2)/eta^2)*(1+tan(LAMDA)^2/Beta^2)))*7.18;

C_L = C_L_alpha*(alpha - alpha_0);

C_D0 = .05; %C_D0 = C_fe*(S_wet/S_ref) -- C_fe = .0035 t12.3

v_stall = sqrt(2*w/(rho_sl*C_L_max*S_wing));

q = 0.5 * rho * v_cruise^2;
q_stall = 0.5 * rho_sl * v_stall^2;

ws_cruise = q*(sqrt((C_D0*pi*AR*e)/3));
ws_loiter = q*(sqrt(pi*AR*e*C_D0));
ws_stall = q_stall*C_L_max;

w_cruise = w*.97*.985;
w_loiter = w_cruise*.856*.9404;
w_stall = w_loiter*.995;

ws_cruise_corr = ws_cruise * (w/w_cruise);
ws_loiter_corr = ws_loiter * (w/w_loiter);
ws_stall_corr = ws_stall * (w/w_stall);

S_wing_new = w/ws_stall_corr;


TW_cruise = (q*C_D0/ws_cruise_corr)+((ws_cruise_corr)/((pi*AR*e)*q));
TW_stall = (q_stall*C_D0/ws_stall_corr)+((ws_stall_corr)/((pi*AR*e)*q_stall));