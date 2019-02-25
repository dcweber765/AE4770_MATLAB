clear variables; close all; clear windows; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
%%%%Inputs%%%%
enginefail = 0; %%0 if no failure, 1 if 1 engine fails
%%FROM XFLR5
dCL_dalpha_deg = 0.0924; %%lift slope of wing and tail
alpha_0L_deg = -5.595; %%Angle of 0 Lift of wing airfoil
airfoil_liftslope_rad = 6.13; %%lift slope of wing airfoil (1/rad) %%CHANGE W/VELOCITY???
airfoil_radome_liftslope = 0.1083 * 180/pi; %%(1/rad)
alpha_0L_radome = -1.801 * pi/180;
%%Plane
w_0 = 40000;
S_w = 627; %%ft^2
b = 76; %%ft
c_root = 12.5; %%ft
mac_wing = 8.98;
fuse_diameter_nose = 6; %%ft
fuse_diameter_max = 8;
fuse_diameter_tail = 6;
fuse_length = 55; %ft
u = 8 * pi/180; %%upsweep angle
sweep_LE_deg = 10.44;
sweep_maxt_deg = 4.9; %%%sweep at chord location where airfoil is thickest
S_ht = 232;
b_ht = 32;
c_root_ht = 10;
mac_ht = 7.6;
S_vt_main = 60;
mac_vt_main = 6.89;
S_vt_side = (33 + 44)/2;
mac_vt_side = 3.12;
radar_diameter = 24;
radar_height = 4;
length_nacelle = 14; %%inches to feet
diameter_nacelle = 3.5; %%inches to feet
diameter_missilestore = 8.5 / 12; %%inches to feet
tc_w = 0.149; %%airfoil thickness/chord
x_cbar_w = 0.44; %%airfoil max thickness location/chord
tc_tail = 0.09;
x_cbar_tail = 0.30;
sweep_LE_ht_deg = 14.04;
sweep_maxt_ht_deg = 8.36;
sweep_LE_vt_deg = 57.99; %%main has max sweep
sweep_maxt_vt_deg = 29.25;
%%Mission Segment Weights
w_cruise = 40000*0.9554; %%lbs
%%Thrust/Power
%%Cruise Flight
alt_cruise = 30000; %%ft
v_cruise_mph = 350;
alpha_deg = 2;
alpha_radome = 0 * pi/180; %%deg to rad
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
%%%%Initial Calcs%%%%
%%Conversions
dCL_dalpha_rad = dCL_dalpha_deg * 180/pi;
alpha_0L_rad = alpha_0L_deg * pi/180;
v_cruise = v_cruise_mph * 1.4666667; %%ft/s
alpha_rad = alpha_deg * pi/180;
sweep_LE = sweep_LE_deg * pi/180;
sweep_maxt_rad = sweep_maxt_deg * pi/180;
sweep_LE_ht = sweep_LE_ht_deg * pi/180;
sweep_maxt_ht = sweep_maxt_ht_deg * pi/180;
sweep_LE_vt = sweep_LE_vt_deg * pi/180;
sweep_maxt_vt = sweep_maxt_vt_deg * pi/180;
%%Plane Areas
AR_w = b^2/S_w;
AR_ht = b_ht^2/S_ht;
S_exposed = S_w - fuse_diameter_max*c_root; %%%IMPROVE THIS!!!!!********
S_exposed_ht = S_ht - fuse_diameter_tail*c_root_ht;
S_radar = 0.25 * pi * radar_diameter^2;
S_wet_wing = 2*S_w;
S_wet_ht = 2*S_ht;
S_vt = S_vt_main + 2 * S_vt_side; %%Total planform
S_wet_vt = 2*S_vt;
S_wet_radar = S_radar * 2;
S_wet_fuse = pi * fuse_diameter_max * fuse_length;
S_wet_nacelle = pi * diameter_nacelle * length_nacelle;
S_wet = S_wet_wing + S_wet_ht + S_wet_vt + S_wet_fuse + 2*S_wet_nacelle + S_wet_radar;
AR_radome = radar_diameter^2/S_radar;
%%Cruise Flight
if alt_cruise > 36090 %%ft, rho in slugs/ft^3
 [temp, rho, mu] = atmosphere_eng_highalt( alt_cruise );
else
 [temp, rho, mu] = atmosphere_eng_lowalt( alt_cruise );
end
q_cruise = 0.5 * rho * v_cruise^2;
M_cruise = v_cruise_mph / 678;
WS_cruise = w_cruise / S_w;
%%Variable
syms V_mph
V_syms = V_mph * 1.4666667;
q_syms = 0.5 * rho * V_syms^2;
M_syms = V_mph / 678;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
%%%%Thrust-Drag Plot%%%%
%%Step 2: CL estimate
CL = dCL_dalpha_rad * (alpha_rad - alpha_0L_rad);
beta_sqrd = 1 - M_cruise^2;
beta_sqrd_syms = 1 - M_syms^2;
nu = airfoil_liftslope_rad / (2*pi/sqrt(beta_sqrd)); %%SHOULD BE NEAR 0.95!!!
nu_syms = airfoil_liftslope_rad / (2*pi/sqrt(beta_sqrd_syms));
%nu = 0.95;
%nu_syms = 0.95;
%F = 1.07 * (1 + fuse_diameter_max/b)^2;
dCL_dalpharad_theoretical = 2*pi*AR_w * 1 / (2 + sqrt(4 + (AR_w^2 * beta_sqrd / nu^2) * ( 1 +(tan(sweep_maxt_rad))^2 / beta_sqrd) ));
dCL_dalpharad_theoretical_syms = 2*pi*AR_w * 1 / (2 + sqrt(4 + (AR_w^2 * beta_sqrd_syms /nu_syms^2) * ( 1 + (tan(sweep_maxt_rad))^2 / beta_sqrd_syms) ));
CL_theoretical = dCL_dalpharad_theoretical * (alpha_rad - alpha_0L_rad);
CL_theoretical_syms = dCL_dalpharad_theoretical_syms * (alpha_rad - alpha_0L_rad);
% %H tailnu_ht = ( 0.113 * 180/pi) / (2*pi/sqrt(beta_sqrd));
nu_ht = 0.95;
%F_ht = 1.07 * (1 + fuse_diameter_tail/b_ht)^2;
dCL_dalpharad_tail = 2*pi*AR_ht * 0.98 / (2 + sqrt(4 + (AR_ht^2 * beta_sqrd / nu_ht^2) * ( 1 +(tan(sweep_maxt_ht))^2 / beta_sqrd) ));
dCL_dalpharad_tail_syms = 2*pi*AR_ht * 0.98 / (2 + sqrt(4 + (AR_ht^2 * beta_sqrd_syms /nu_ht^2) * ( 1 + (tan(sweep_maxt_ht))^2 / beta_sqrd_syms) ));
%Radome
nu_radome = ( airfoil_radome_liftslope) / (2*pi/sqrt(beta_sqrd));
nu_radome_syms = ( airfoil_radome_liftslope) / (2*pi/sqrt(beta_sqrd_syms));
%nu_radome = 0.95;
%nu_radome_syms = 0.95;
%F_radome = 1.07 * (1 + 1.5/radar_diameter)^2;
dCL_dalpharad_radome = 2*pi*AR_radome * 1 / (2 + sqrt(4 + (AR_radome^2 * beta_sqrd /nu_radome^2) * ( 1 + (tan(0))^2 / beta_sqrd) ));
dCL_dalpharad_radome_syms = 2*pi*AR_radome * 1 / (2 + sqrt(4 + (AR_radome^2 *beta_sqrd_syms / nu_radome_syms^2) * ( 1 + (tan(0))^2 / beta_sqrd_syms) ));
CL_radome = dCL_dalpharad_radome * (alpha_radome - alpha_0L_radome);
CL_radome_syms = dCL_dalpharad_radome_syms * (alpha_radome - alpha_0L_radome);
%Est. Total
% dCL_dalpharad_total = dCL_dalpharad_radome*S_radar/S_w + dCL_dalpharad_theoretical;
% dCL_dalpharad_total_syms = ( dCL_dalpharad_radome_syms*S_radar +dCL_dalpharad_theoretical_syms*S_w ) / (S_radar + S_w);
%
% CL_weighted_total = ( CL_radome*S_radar + CL_theoretical*S_w ) / (S_radar + S_w);
% CL_weighted_total_syms = ( CL_radome_syms*S_radar + CL_theoretical_syms*S_w ) /
(S_w);
%
% CL_equivalent = CL_weighted_total * (S_radar + S_w)/S_w;
% CL_equivalent_syms = CL_weighted_total_syms * (S_radar + S_w)/S_w;
% CL_equivalent = CL_radome*S_radar/S_w + CL_theoretical;
% CL_equivalent_syms = CL_radome_syms*S_radar/S_w + CL_theoretical_syms;
%%Total
taper_w = c_root/4;
sweep_qc = 7.31 * pi/180;
l_bar_vt = 22.4513;
l_bar_ht = 22.6433;
de_dalpha = 4.44 * sqrt(1-M_cruise^2) * ( ( 1/AR_w - 1/(1+AR_w^1.7) ) * (10 - 3*taper_w)/7 * (1- l_bar_vt/b)/(2*l_bar_ht/b)^0.33 * sqrt(cos(sweep_qc)) )^1.19;
de_dalpha_syms = 4.44 * sqrt(1-M_syms^2) * ( ( 1/AR_w - 1/(1+AR_w^1.7) ) * (10 -3*taper_w)/7 * (1 - l_bar_vt/b)/(2*l_bar_ht/b)^0.33 * sqrt(cos(sweep_qc)) )^1.19;
a_bar = dCL_dalpharad_theoretical * ( 1 +(dCL_dalpharad_radome/dCL_dalpharad_theoretical) * (S_radar/S_w) );
a_bar_syms = dCL_dalpharad_theoretical_syms * ( 1 +(dCL_dalpharad_radome_syms/dCL_dalpharad_theoretical_syms) * (S_radar/S_w) );
CL_equivalent = a_bar * (alpha_rad - alpha_0L_rad);
CL_equivalent_syms = a_bar_syms * (alpha_rad - alpha_0L_rad);
L_cruise = q_cruise * S_w * CL_equivalent;
L_cruise_syms = q_syms * S_w * CL_equivalent_syms;
%%Step 3: CD_0 (Parastite/Zero-Lift Drag) Estimate
CD_0_equivSF = 0.0035 * S_wet/S_w;
CD_0_equivSF2 = 0.003 * S_wet/S_w;
%Components
R_wing = rho * v_cruise * mac_wing / mu;
R_wing_syms = rho * V_syms * mac_wing / mu;
Cf_wing = 0.455 / ( (log10(R_wing))^2.58 * (1 + 0.144*M_cruise^2)^0.65 );
Cf_wing_syms = 0.455 / ( (log10(R_wing_syms))^2.58 * (1 + 0.144*M_syms^2)^0.65 );
FF_wing = ( 1+( 0.6*tc_w/x_cbar_w )+( 100*(tc_w^4) ) ) * ( 1.34*(M_cruise^0.18) *(cos(sweep_maxt_rad))^0.28 ); %Form Factor Shape
FF_wing_syms = ( 1+( 0.6*tc_w/x_cbar_w )+( 100*(tc_w^4) ) ) * ( 1.34*(M_syms^0.18) *(cos(sweep_maxt_rad))^0.28 ); %Form Factor Shape
CBU_wing = Cf_wing * FF_wing * S_wet_wing;
CBU_wing_syms = Cf_wing_syms * FF_wing_syms * S_wet_wing;
R_ht = rho * v_cruise * mac_ht / mu;
R_ht_syms = rho * V_syms * mac_ht / mu;
Cf_ht = 0.455 / ( (log10(R_ht))^2.58 * (1 + 0.144*M_cruise^2)^0.65 );
Cf_ht_syms = 0.455 / ( (log10(R_ht_syms))^2.58 * (1 + 0.144*M_syms^2)^0.65 );
FF_ht = ( 1+( 0.6*tc_tail/x_cbar_tail )+( 100*(tc_tail^4) ) ) * ( 1.34*(M_cruise^0.18) *(cos(sweep_maxt_ht))^0.28 ); %Form Factor Shape
FF_ht_syms = ( 1+( 0.6*tc_tail/x_cbar_tail )+( 100*(tc_tail^4) ) ) * ( 1.34*(M_syms^0.18) *(cos(sweep_maxt_ht))^0.28 ); %Form Factor Shape
Q_ht = 1.08;
CBU_ht = Cf_ht * FF_ht * S_wet_ht * Q_ht;
CBU_ht_syms = Cf_ht_syms * FF_ht_syms * S_wet_ht * Q_ht;
R_vt_main = rho * v_cruise * mac_vt_main / mu;
R_vt_main_syms = rho * V_syms * mac_vt_main / mu;
Cf_vt_main = 0.455 / ( (log10(R_vt_main))^2.58 * (1 + 0.144*M_cruise^2)^0.65 );
Cf_vt_main_syms = 0.455 / ( (log10(R_vt_main_syms))^2.58 * (1 + 0.144*M_syms^2)^0.65 );
FF_vt_main = ( 1+( 0.6*tc_tail/x_cbar_tail )+( 100*(tc_tail^4) ) ) * ( 1.34*(M_cruise^0.18) *(cos(sweep_maxt_vt))^0.28 ); %Form Factor Shape
FF_vt_main_syms = ( 1+( 0.6*tc_tail/x_cbar_tail )+( 100*(tc_tail^4) ) ) * ( 1.34*(M_syms^0.18) *(cos(sweep_maxt_vt))^0.28 ); %Form Factor Shape
Q_vt = 1.08;
CBU_vt_main = Cf_vt_main * FF_vt_main * 2*S_vt_main * Q_vt;
CBU_vt_main_syms = Cf_vt_main_syms * FF_vt_main_syms * 2*S_vt_main * Q_vt;
R_vt_side = rho * v_cruise * mac_vt_side / mu;
R_vt_side_syms = rho * V_syms * mac_vt_side / mu;
Cf_vt_side = 0.455 / ( (log10(R_vt_side))^2.58 * (1 + 0.144*M_cruise^2)^0.65 );
Cf_vt_side_syms = 0.455 / ( (log10(R_vt_side_syms))^2.58 * (1 + 0.144*M_syms^2)^0.65 );
FF_vt_abv = ( 1+( 0.6*tc_tail/x_cbar_tail )+( 100*(tc_tail^4) ) ) * ( 1.34*(M_cruise^0.18) *(cos(17.03*pi/180))^0.28 ); %Form Factor Shape
FF_vt_abv_syms = ( 1+( 0.6*tc_tail/x_cbar_tail )+( 100*(tc_tail^4) ) ) * ( 1.34*(M_syms^0.18) *(cos(17.03*pi/180))^0.28 ); %Form Factor Shape
FF_vt_blw = ( 1+( 0.6*tc_tail/x_cbar_tail )+( 100*(tc_tail^4) ) ) * ( 1.34*(M_cruise^0.18) *(cos(22.21*pi/180))^0.28 ); %Form Factor Shape
FF_vt_blw_syms = ( 1+( 0.6*tc_tail/x_cbar_tail )+( 100*(tc_tail^4) ) ) * ( 1.34*(M_syms^0.18) *(cos(22.21*pi/180))^0.28 ); %Form Factor Shape
Q_vt = 1.08;
CBU_vt_abv = Cf_vt_side * FF_vt_abv * 2*44 * Q_vt;
CBU_vt_abv_syms = Cf_vt_side_syms * FF_vt_abv_syms * 2*44 * Q_vt;
CBU_vt_blw = Cf_vt_side * FF_vt_blw * 2*33 * Q_vt;
CBU_vt_blw_syms = Cf_vt_side_syms * FF_vt_blw_syms * 2*33 * Q_vt;
R_fuse = rho * v_cruise * fuse_length / mu;
R_fuse_syms = rho * V_syms * fuse_length / mu;
Cf_fuse = 0.455 / ( (log10(R_fuse))^2.58 * (1 + 0.144*M_cruise^2)^0.65 );
Cf_fuse_syms = 0.455 / ( (log10(R_fuse_syms))^2.58 * (1 + 0.144*M_syms^2)^0.65 );
f_fuse = fuse_length/fuse_diameter_max;
FF_fuse = 1 + 60/f_fuse^3 + f_fuse/400;
CBU_fuse = Cf_fuse * FF_fuse * S_wet_fuse;
CBU_fuse_syms = Cf_fuse_syms * FF_fuse * S_wet_fuse;
R_nacelle = rho * v_cruise * length_nacelle / mu;
R_nacelle_syms = rho * V_syms * length_nacelle / mu;
Cf_nacelle = 0.455 / ( (log10(R_nacelle))^2.58 * (1 + 0.144*M_cruise^2)^0.65 );
Cf_nacelle_syms = 0.455 / ( (log10(R_nacelle_syms))^2.58 * (1 + 0.144*M_syms^2)^0.65 );
f_nacelle = length_nacelle/diameter_nacelle;
FF_nacelle = 1 + 0.35/f_nacelle;
Q_nacelle = 1.5;
CBU_nacelle = Cf_nacelle * FF_nacelle * S_wet_nacelle * Q_nacelle;
CBU_nacelle_syms = Cf_nacelle_syms * FF_nacelle * S_wet_nacelle * Q_nacelle;
R_radar = rho * v_cruise * radar_diameter / mu;
R_radar_syms = rho * V_syms * radar_diameter / mu;
Cf_radar = 0.455 / ( (log10(R_radar))^2.58 * (1 + 0.144*M_cruise^2)^0.65 );
Cf_radar_syms = 0.455 / ( (log10(R_radar_syms))^2.58 * (1 + 0.144*M_syms^2)^0.65 );
f_radar = 1;
FF_radar = 1 + 0.35/f_radar;
CBU_radar = Cf_radar * FF_radar * S_wet_radar;
CBU_radar_syms = Cf_radar_syms * FF_radar * S_wet_radar;
%Misc. Drag
Dq_missilestore = 0.15;
FA_missilestore = 0.5 * 0.25 * pi * diameter_missilestore^2;
CD_missilestore = Dq_missilestore / S_w;
Dq_upsweep = 3.83 * ( 0.25 * pi * fuse_diameter_max^2 ) * u^2.5;
CD_upsweep = Dq_upsweep / S_w;
Dq_radarsupport = 0.3;
CD_radarsupport = Dq_radarsupport / S_w;
Dq_strut = 0.05;
CD_struts = Dq_strut / S_w;
if enginefail == 0
 CD_engfail = 0;
else 
 CD_engfail = 0.1 * (0.04*6) * 0.25 * pi * 14^2 / S_w;
end
CD_misc = CD_missilestore + CD_engfail + CD_upsweep + CD_struts + CD_radarsupport;
CD_LP_fraction = 1.05;
%Totals
CD_02 = ( ( (CBU_wing + CBU_ht + CBU_vt_main + 2*CBU_vt_abv + 2*CBU_vt_blw +CBU_fuse + 2*CBU_nacelle + CBU_radar)/S_w ) + CD_misc ) * CD_LP_fraction;
CD_02_syms = ( ( (CBU_wing_syms + CBU_ht_syms + CBU_vt_main_syms +2*CBU_vt_abv_syms + 2*CBU_vt_blw_syms + CBU_fuse_syms + 2*CBU_nacelle_syms +CBU_radar_syms)/S_w ) + CD_misc ) * CD_LP_fraction;
CD_0 = ( ( (CBU_wing + CBU_ht + CBU_vt_main + 2*CBU_vt_abv + 2*CBU_vt_blw +CBU_fuse + 2*CBU_nacelle + CBU_radar)/S_w ) );
CD_0_syms = ( ( (CBU_wing_syms + CBU_ht_syms + CBU_vt_main_syms +2*CBU_vt_abv_syms + 2*CBU_vt_blw_syms + CBU_fuse_syms + 2*CBU_nacelle_syms +CBU_radar_syms)/S_w ) );
%%Step 4: CD_ind (Induced Drag) Estimate
e_0_30degsweep = 4.61 * (1 - 0.045*AR_w^0.68) * (cos(30*pi/180))^0.15 - 3.1;
e_0_0sweep = 1.78 * (1 - 0.045*AR_w^0.68) - 0.64;
e_0 = (e_0_30degsweep - e_0_0sweep) * sweep_LE_deg / 30 + e_0_0sweep;
% CD_indw_est = CL^2 / ( pi * AR_w * e_0 );
% CD_ind_theoreticalw = CL_theoretical^2 / ( pi * AR_w * e_0 );
CD_ind_theoretical = CL_equivalent^2 / ( pi * AR_w * e_0 );
CD_ind_theoretical_syms = CL_equivalent_syms^2 / ( pi * AR_w * e_0 );
% CD_ind_radome = (CL_radome * S_radar/S_w)^2 / ( pi * 1 * e_0 );
% CD_ind_equivalent = CD_ind_theoreticalw + CD_ind_radome;
mu = radar_diameter/b;
r = (CL_radome * S_radar) / (CL_theoretical * S_w);
r_syms = (CL_radome_syms * S_radar) / (CL_theoretical_syms * S_w);
gap_avgspan = radar_height / ( 0.5 * (radar_diameter + b) );
sigma = 0.1;
e_0_biplane = 0.8 * mu^2 * (1 + r)^2 / (mu^2 + 2 * sigma * mu * r + r^2);
e_0_biplane_syms = 0.8 * mu^2 * (1 + r_syms)^2 / (mu^2 + 2 * sigma * mu * r_syms +r_syms^2);
CD_ind_biplane = CL_equivalent^2 / ( pi * AR_w * e_0_biplane );
CD_ind_biplane_syms = CL_equivalent_syms^2 / ( pi * AR_w * e_0_biplane );
% CD_ind = CD_ind_theoretical;
% CD_ind_syms = CD_ind_theoretical_syms;
%Assume L = W
D_ind = w_cruise^2/(e_0 * pi * AR_w * S_w * q_cruise);
D_ind_syms = w_cruise^2/(e_0 * pi * AR_w * S_w * q_syms);
CD_ind = D_ind/(q_cruise * S_w);
%%Step 6: Drag
D_0 = (CD_0 + CD_engfail) * q_cruise * S_w;
D_0_syms = CD_0_syms * q_syms * S_w;
D = D_0 + D_ind;
D_syms = D_0_syms + D_ind_syms;
% figure()
% fplot([D_syms, D_0_syms, D_ind_syms], [120 450])
% title("Drag Plot")
% xlabel("Velocity (mph)")
% ylabel("Drag Force (lbs)")
% legend("Total Drag", "Parasitic Drag", "Induced Drag")
%%Step 7: Thrust
TW_cruise = q_cruise * CD_02 / WS_cruise + WS_cruise / ( pi * AR_w * e_0 * q_cruise );
TW_syms = q_syms * CD_02_syms / WS_cruise + WS_cruise / ( pi * AR_w * e_0 * q_syms );
T_cruise = TW_cruise * w_cruise;
T_syms = TW_syms * w_cruise;
T_max = 2 * 550 * 5071 * 0.9 / v_cruise;
T_max_syms = 2 * 550 * 5071 * 0.9 / V_syms;
WScruise_maxR = q_cruise * sqrt( pi * AR_w * e_0 * CD_0 );
WScruise_maxR_syms = q_syms * sqrt( pi * AR_w * e_0 * CD_0 );
T_maxR = ( q_cruise * CD_0 / WScruise_maxR + WScruise_maxR / ( pi * AR_w * e_0 *q_cruise ) ) * w_cruise;
T_maxR_syms = ( q_syms * CD_0 / WScruise_maxR_syms + WScruise_maxR_syms / ( pi *AR_w * e_0 * q_syms ) ) * w_cruise;
% figure()
% fplot([T_syms], [100 450])
% title("Thrust-Drag Plot")
% xlabel("Velocity (mph)")
% ylabel("Thrust and Drag Forces (lbs)")
%%RESULTS
TD_cruise = T_cruise/D;
L_shouldbe = w_cruise;
LD_shouldbe = L_shouldbe/D;
LD = L_cruise/D;
TW_inv = 1/TW_cruise;
% figure()
% fplot([T_cruise, D_syms], [100 425])
% title("Thrust-Drag Plot")
% xlabel("Velocity (mph)")
% ylabel("Thrust and Drag Forces (lbs)")
% legend("Req. Thrust for Cruise at 350 mph", "Drag")
%
% figure()
% fplot([T_cruise, T_maxR, D_syms], [100 425])
% title("Thrust-Drag Plot")
% xlabel("Velocity (mph)")
% ylabel("Thrust and Drag Forces (lbs)")
% legend("Req. Thrust for Cruise at 350 mph", "Req. Thrust for Max. Range", "Drag")