
clc;
clear all

mach = 1.5;
w = 1.6028e+04;
%b = 30;%ft span
S_wet = 354;
S_wing = w/50;


AR_eqv = 4.737*(mach^-.979);
lamda = .3;
b = sqrt(AR_eqv*S_wing);
c_root = (2*S_wing)/(b*(1+lamda));
c_tip = lamda*c_root;
wing_sweep = 40;
%AR = (b^2)/(S_wing)



rho = 8.899028e-04;%@30000 ft

v_crusie = 1125*.9;

q = .5*rho*v_crusie^2;

Coe_ht = .7;
Coe_vt = .06;

L_ht = 15;
L_vt = 12;

lamda = .3;
b = sqrt(AR_eqv*S_wing);
c_root = (2*S_wing)/(b*(1+lamda));
c_tip = lamda*c_root;

c_bar = (2/3)*c_root*(1+lamda+lamda^2)/(1+lamda);
Y_bar = (b/6)*((1+2*lamda)/(1+lamda));

S_ht = (Coe_ht*S_wing*c_bar)/L_ht;
S_vt = (Coe_vt*S_wing*b)/L_vt;

CL_req = (w/(.5*rho*v_crusie^2*S_wing))*(1+(2/AR_eqv));
CL = w * (1 + 2/AR_eqv) / (q * S_wing)

a_fuse = .79;%t6.3
Coe_fuse = .41;

fuse_length = a_fuse*w^Coe_fuse;
fineness = 14;%pg 157
fuse_thick = fineness/fuse_length;

AR_ht = 3.5;

b_ht = sqrt(AR_ht*S_ht);
c_root_ht = (2*S_ht)/(b_ht*(1+lamda));
c_tip_ht = lamda*c_root_ht;
c_bar_ht = (2/3)*c_root_ht*(1+lamda+lamda^2)/(1+lamda);
Y_bar_ht = (b_ht/6)*((1+2*lamda)/(1+lamda));

AR_vt = 1/.7;

b_vt = sqrt(AR_vt*S_vt);
c_root_vt = (2*S_vt)/(b_vt*(1+lamda));
c_tip_vt = lamda*c_root_vt;
c_bar_vt = (2/3)*c_root_vt*(1+lamda+lamda^2)/(1+lamda);
Y_bar_vt = (b_vt/6)*((1+2*lamda)/(1+lamda));


