clear all;

MACH = 1.5;
S_ref = 320.56;
S_wing = 232.6;
AR = 4.737*(MACH^-.979);
h = 30000;
LAMDA_LE = 40;
w= 16028;

g = 32.15;                          %Gravitational acceleration [ft/sec^2]
R = 3089.2;                         %Specific gas constant for dry air [lb*ft/sl*K]
a1 = -1.9812*10^-3;                 %Temperature gradient [K/ft]
rho_sl = 0.0023769;                 %Rho at sea level [sl/ft^3]
tau_sl = 288.16;                    %Standard temperature at sea level [K]
mu_sl = 3.737E-7;                   %Dynamic Viscostiy at sea level [lb s/ft2]
S_mu = 110.4;                       %Sutherland Constant for Viscosity [K}



if h<=36089.2                       %For altitudes 0 - 36089.2 km
    tau = tau_sl + a1*h;
    rho = rho_sl * (tau/tau_sl).^(-1-g/(a1*R));
    mu = mu_sl*(tau/tau_sl)^3/2 * ((tau_sl + S_mu)/ (tau + S_mu));
elseif h>36089.2                    %For altitude 36089.2 - 82021 km
    tau_11 = tau_sl + a1*36089.2;   %Temperature at 36089.2 km
    rho_11 = rho_sl * (tau/tau_sl).^(-1-g/(a1*R));    %Pressure at 36089.2 km
    
    tau = tau_11;
    rho = rho_11*exp((-g*(h-36089.2))/(R*tau_11));
    mu = mu_sl*(tau/tau_sl)^3/2 * ((tau_sl + S_mu)/ (tau + S_mu));
end


k = 2.08E-5;%Smooth Paint Drag



v_range = linspace(100,1100);

for i = 1:size(v_range,2)
i
v  = v_range(i)
mach = v/1125
%% Induced Drag

q = .5*rho*v^2;

e = 4.61*(1-0.045*AR^.68)*(cosd(LAMDA_LE))^.15 - 3.1;

K = 1/(pi*AR*e);

K_super = (AR*(mach^2 - 1)*cosd(40))/(4*AR*sqrt(mach^2 - 1)) - 2;


D_Induced(i) = (w^2) / (q*pi*AR*S_wing*e);%K*(w^2)/q;


%% Paracitic Drag

%wings
l_wings = 11;
Re_wings = (rho*v*l_wings)/mu;
R_cutoff_wings = 44.62*(l_wings/k)^1.053 * mach^1.16;
% if(Re_wings<R_cutoff_wings)
%     Re_wings = R_cutoff_wings;
% end

C_f_wings = 0.455/(log10(Re_wings)^2.58 * (1 + .144*mach^2)^.65);
%MS313
xc_ms313 = .375;
tc_ms313 = .131;
FF_wings = (1 + (.6/xc_ms313)*tc_ms313 + 100*tc_ms313^4)*(1.34*mach^.18*(cosd(40))^.28);
S_wet_wings = 233.6391;


%HT
l_HT = 7.52;
Re_HT = (rho*v*l_HT)/mu;
R_cutoff_HT = 44.62*(l_HT/k)^1.053 * mach^1.16;
% if(Re_HT<R_cutoff_HT)
%     Re_HT = R_cutoff_HT;
% end
C_f_HT = 0.455/(log10(Re_wings)^2.58 * (1 + .144*mach^2)^.65);
%0012
xc_0012 = .3;
tc_0012 = .12;
FF_HT = (1 + (.6/xc_0012)*tc_0012 + 100*tc_0012^4)*(1.34*mach^.18*(cosd(45))^.28);
S_wet_HT = 124.4474;

%VT
l_VT = 6.566;
Re_VT = (rho*v*l_VT)/mu;
R_cutoff_VT = 44.62*(l_VT/k)^1.053 * mach^1.16;
% if(Re_VT<R_cutoff_VT)
%     Re_VT = R_cutoff_VT;
% end
C_f_VT = 0.455/(log10(Re_wings)^2.58 * (1 + .144*mach^2)^.65);
%0012
xc_0012 = .3;
tc_0012 = .12;
FF_VT = (1 + (.6/xc_0012)*tc_0012 + 100*tc_0012^4)*(1.34*mach^.18*(cosd(40))^.28);
S_wet_VT = 63.6083;

%Fusealge
l_fuse = 41.8;
A_max_fuse = 25;
Re_fuse = (rho*v*l_fuse)/mu;
R_cutoff_fuse = 44.62*(l_fuse/k)^1.053 * mach^1.16;
% if(Re_fuse<R_cutoff_fuse)
%     Re_fuse = R_cutoff_fuse;
% end
C_f_fuse = 0.455/(log10(Re_fuse)^2.58 * (1 + .144*mach^2)^.65);
f_fuse = l_fuse/(sqrt((4/pi)*A_max_fuse));
FF_fuse = (.9 + (5/f_fuse^1.5) + (f_fuse/400));
S_wet_fuse = 735.1544;


%Canopy
l_canopy = 12;
A_max_canopy = 5.5;
Re_canopy = (rho*v*l_canopy)/mu;
R_cutoff_canopy = 44.62*(l_canopy/k)^1.053 * mach^1.16;
% if(Re_canopy<R_cutoff_canopy)
%     Re_canopy = R_cutoff_canopy;
% end
C_f_canopy = 0.455/(log10(Re_canopy)^2.58 * (1 + .144*mach^2)^.65);
f_canopy = l_canopy/(sqrt((4/pi)*A_max_canopy));
FF_canopy = (.9 + (5/f_canopy^1.5) + (f_canopy/400));
S_wet_canopy = 66.2804;


%Inlet
l_inlet = 3;
A_max_inlet = pi*2.1355^2;
Re_inlet = (rho*v*l_inlet)/mu;
R_cutoff_inlet = 44.62*(l_inlet/k)^1.053 * mach^1.16;
% if(Re_inlet<R_cutoff_inlet)
%     Re_inlet = R_cutoff_inlet;
% end
C_f_inlet = 0.455/(log10(Re_inlet)^2.58 * (1 + .144*mach^2)^.65);
f_inlet = l_inlet/(sqrt((4/pi)*A_max_inlet));
FF_inlet = 1 + (.35/f_inlet);
S_wet_inlet = 21.1441;
Q_inlet = 1.5;


% Wave Drag
A_max_total = 43.2;
lentgh_term = 46;
Dq_SH = (9*pi/2)*(A_max_total/lentgh_term)^2;

E_wd = 2.1;

Dq_wd = E_wd*(1-.386*(mach-1.2)^.57 * (1-((pi*40^.77) / 100)))*Dq_SH;

if (v > 1125)
    FF_wings = 1;
    FF_HT= 1;
    FF_VT= 1;
    FF_fuse= 1;
    FF_canopy= 1;
    FF_inlet= 1;
    C_D_wave = Dq_wd;
else
    C_D_wave = 0;
end

C_D0_wings = C_f_wings*FF_wings*S_wet_wings; 
C_D0_HT = C_f_HT*FF_HT*S_wet_HT;
C_D0_VT = C_f_VT*FF_VT*S_wet_VT;
C_D0_fuse = C_f_fuse*FF_fuse*S_wet_fuse;
C_D0_canopy = C_f_canopy*FF_canopy*S_wet_canopy;
C_D0_inlet = C_f_inlet*FF_inlet*S_wet_inlet*Q_inlet;

S_wet = 2*S_wet_wings + 2*S_wet_HT + 2*S_wet_VT + 2*S_wet_fuse + S_wet_canopy + 2*S_wet_inlet;

C_D_misc = 0;
C_D_LP = 0;

C_D0 =  ((2*C_D0_wings + C_D0_fuse + C_D0_canopy + 2*C_D0_inlet + 2*C_D0_HT +2*C_D0_VT)/S_wing) + C_D_misc;

D0_i(i) = q*C_D0*S_wet;
end

D0_total= D0_i + D_Induced;

plot(v_range, D0_i, 'red')
hold on
plot(v_range, D_Induced, 'blue')
plot(v_range, D0_total ,'green')
hold off
title('Total Drag')
xlabel('Velocity [ft/s]')
ylabel('Drag Force [lbs]')

figure
plot(v_range, D0_total)
hold on
thrust = 21000*ones(1,100);
plot(v_range, thrust) 
