
clc;
clear all

mach = 1.5;
b = 30;%ft span
S_wet = 354;
S_wing = 16000/50;
AR = (b^2)/(S_wing)

AR_eqv = 4.737*(mach^-.979)

w = 1.6028e+04;
rho = 8.899028e-04;%@30000 ft

v_crusie = 1125*.9;

q = .5*rho*v_crusie^2;

Coe_ht = .7;
Coe_vt = .06;

L_ht = 12;
L_vt = 9;

c_tip = 4;
c_root = 10;

lamda = c_tip/c_root;

c_bar = (2/3)*c_root*(1+lamda+lamda^2)/(1+lamda);


S_ht = (Coe_ht*S_wing*c_bar)/L_ht;
S_vt = (Coe_vt*S_wing*b)/L_vt;

CL_req = (w/(.5*rho*v_crusie^2*S_wing))*(1+(2/AR_eqv));
CL = w * (1 + 2/AR_eqv) / (q * S_wing)