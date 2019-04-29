clear all;

g = 32.15;                          %Gravitational acceleration [ft/sec^2]
R = 3089.2;                         %Specific gas constant for dry air [lb*ft/sl*K]
a1 = -1.9812*10^-3;                 %Temperature gradient [K/ft]
rho_sl = 0.0023769;                 %Density at sea level [sl/ft^3]
tau_sl = 288.16;                    %Standard temperature at sea level [K]


altt = 30000;               %Altitude in feet

    h = altt;
    
    if h<=36089.2                           %For altitudes 0 - 36089.2 km
        tau = tau_sl + a1*h;
        rho = rho_sl * (tau/tau_sl).^(-1-g/(a1*R));
    elseif h>36089.2                        %For altitude 36089.2 - 82021 km
        tau_11 = tau_sl + a1*36089.2;       %Temperature at 36089.2 km
        rho_11 = rho_sl * (tau_11/tau_sl).^(-1-g/(a1*R));    %Pressure at 36089.2 km

        tau = tau_11;
        rho = rho_11*exp((-g*(h-36089.2))/(R*tau_11));
    end

S = 320.5;                            %Total Wing surface area [ft^2]
S_t = 51.2;
Span = 31.9;                        %Wing Span [ft]
CBar = 11;                       %Mean Aerodynamic wing chord
C_D0 = 0.02;                       %Zero-lift drag coefficient [unitless]
W = 1.6028e+04;                          %Weight of the aircraft [lb]
mach = 1.5;
AR = 4.737*(mach^-.979);                          %Aspect ratio [unitless]
LAMDA_LE = 40;
e = 4.61*(1-0.045*AR^.68)*(cosd(LAMDA_LE))^.15 - 3.1;                            %Oswald efficience factor [unitless]
h_ac = .25;
a = 6.06;
a_t = 4.7;
a_e = .4;
i_t = deg2rad(6);
dEpdAl = .17;
l_t = 15;
V_H = (S_t/S)*(l_t/CBar);
aBar = a + a_t*(1-dEpdAl)*(S_t/S);
V = 1012.5;
%rho =  0.0018;
Q = .5*rho*(V^2);
epsilon_0 = 0;
C_l0w = .1;
C_mac = .01;
h = .1*CBar;



h_n = h_ac + (a_t/aBar)*V_H*(1-dEpdAl);

C_l0 = C_l0w + (S_t/S)*a_t*(i_t - epsilon_0);
C_m0 = C_mac + C_l0w*(h-h_ac) - a_t*(i_t - epsilon_0)*V_H*(1-(h-h_ac)*(CBar/l_t));


C_lAlpha = aBar;
C_lDeltaE = (S_t/S)*a_e;
C_mAlpha = aBar*(h-h_n) - a_t*V_H*dEpdAl;
C_mDeltaE = C_lDeltaE*(h-h_n) - a_e*V_H;

ClCm = [C_lAlpha C_lDeltaE;...
        C_mAlpha C_mDeltaE];

trim = inv(ClCm)*[((W/(Q*S)) - C_l0); -C_m0]

deg_trim = rad2deg(trim)
