clear variables; close all; clear windows; clc;
%%%Inputs
alpha_deg = 2; %%deg
alt_cruise = 30000; %%ft
w_0 = 160000; %%lbs
v_cruise_mph = .9*746; %%mph
b = 30; %%ft
S_w = 350;%720; %%ft^2
%%%Calculations
if alt_cruise > 36090 %%ft, rho in slugs/ft^3
 [temp, rho] = atmosphere_eng_highalt( alt_cruise );
else
 [temp, rho] = atmosphere_eng_lowalt( alt_cruise );
end

w_cruise = w_0;
v_cruise = v_cruise_mph * 1.4666667; %%ft/s
q = 0.5 * rho * v_cruise^2; %slugs/ft-s^2
AR_w = (b^2)/S_w;
CL_req = w_cruise * (1 + 2/AR_w) / (q * S_w);
Cl_airfoil_req = CL_req;
wingloading_cruise = w_cruise/S_w;
%%%Results
fprintf("Required Cl: " + Cl_airfoil_req + " at " + alpha_deg + " deg \n")
fprintf("Aspect Ratio: " + AR_w + "\n")
fprintf("Cruise Wing Loading: " + wingloading_cruise + " lbs/ft^2 \n\n")

%%%%Functions%%%%
function [ temp_h_eng, rho_h_eng ] = atmosphere_eng_lowalt( alt_ft )
rho_sl = 1.225; %kg/m^3
temp_sl = 288.16; %K
a_1 = -6.5*10^(-3); %K/m
R = 287;
g = 9.8;
alt_m = alt_ft*.3048; %%m
temp_h = temp_sl + a_1*alt_m;
rho_h = rho_sl*(temp_h/temp_sl)^(-1-g/(a_1*R));
temp_h_eng = temp_h;
rho_h_eng = rho_h*0.00194032; %%slugs/ft^3
end
function [ temp_h_eng, rho_h_eng ] = atmosphere_eng_highalt( alt_ft )
rho_sl = 1.225; %kg/m^3
temp_sl = 288.16; %K
a_1 = -6.5*10^(-3); %K/m
R = 287;
g = 9.8;
alt_m = alt_ft*.3048; %%m
temp_11 = temp_sl + a_1*11000;
rho_11 = rho_sl*(temp_11/temp_sl)^(-1-g/(a_1*R));
temp_h = 216.6;
rho_h = rho_11*exp(-g*(alt_m-11000)/(R*temp_11));
temp_h_eng = temp_h;
rho_h_eng = rho_h*0.00194032; %%slugs/ft^3
end