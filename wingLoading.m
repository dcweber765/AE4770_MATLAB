clear all;


%q = 0.5 * rho * v_cruise^2;
ws_crusie = q*(sqrt((C_D0*pi*AR*e)/3));
ws_loiter = q*(sqrt(pi*AR*e*C_D0));

tw_to = .4;
C_L_to = 1.8;

TOP = 150;%fig 5.4
ws_to = C_L_to * tw_to * TOP;