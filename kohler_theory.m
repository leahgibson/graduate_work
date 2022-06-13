% Kohler Theory: Classical and Kappa

%%% The goal of this code is to plot the Kohler curves for two differnet
%%% particles, find the critical supersaturation for each, and find the
%%% Kappa value for each 

clc;
close all;
clear all;

%% parameters 


%%% general parameters %%%

surface_tension_water = 0.073; % J/m^2
Mw = 0.018; % molecular weight of wtaer, kg/mol
density_water = 1000; % kg/m^3
T = 298; % temperature, K
R = 8.314; % universal gas constant, J/molK


%%% particle A: ammonium sulfate %%%

dp_nh4so4 = 50; % dry diameter of ammonium sulfate, nm
density_nh4so4 = 1.77; % density of ammonium sulfate, g/cm^3
Mnh4so4 = 0.132; % molar mass ammonium sulfate, kg/mol
v_nh4so4 = 3; % van't Hoff factor


%%% particle B: sodium chloride %%%

dp_nacl = 50; % dry diameter of sodium chloride, nm
density_nacl = 2.16; % density of sodium chloride, g/cm^3
Mnacl = 0.058; % molar mass sodium chloride, kg/mol
v_nacl = 2; % van't Hoff factor


%%% conversions, SI units %%%

dp_nh4so4 = dp_nh4so4/(1e9); % m
density_nh4so4 = density_nh4so4*((1e2)^3/1000); % kg/m^3
dp_nacl = dp_nacl/(1e9); % m 
density_nacl = density_nacl*((1e2)^3/1000); % kg/m^3

%% Setup mesh

dp = dp_nacl;
n = 1000;
Dp = linspace(dp,1e-5,n); % wet diameters



%% Kohler Curve for particle A


A_nh4so4 = (4 * Mw * surface_tension_water)/(R * T * density_water);
moles_nh4so4 = (v_nh4so4 * pi * dp_nh4so4^3 * density_nh4so4)/(6*Mnh4so4);
B_nh4so4 = (6 * moles_nh4so4 * Mw)/(pi * density_water);


%%% Kolher Curve %%%
Seq_nh4so4 = zeros([1 n]);


for i = 1:n
    Seq_nh4so4(i) = exp(A_nh4so4/Dp(i))*exp(-B_nh4so4/(Dp(i)^3));
end

%%% Critical superaturation %%%
Scrit_nh4so4 = exp(sqrt((4*A_nh4so4^3)/(27*B_nh4so4)))

[M_nh4so4, I_nh4so4] = max(Seq_nh4so4);

%%% Kappa value %%%
Kappa_na4so4 = (4 * A_nh4so4^3)/(27 * dp_nh4so4^3 * (log(Scrit_nh4so4))^2)


%% Kohler Curve for particle B


A_nacl = (4 * Mw * surface_tension_water)/(R * T * density_water);
moles_nacl = (v_nacl * pi * dp_nacl^3 * density_nacl)/(6 * Mnacl);
B_nacl = (6 * moles_nacl * Mw)/(pi * density_water);


%%% Kohler CCurve %%%
Seq_nacl = zeros([1 n]);


for i = 1:n
    Seq_nacl(i) = exp(A_nacl/Dp(i))*exp(-B_nacl/(Dp(i)^3));
end


%%% Critical supersaturation %%%
Scrit_nacl = exp(sqrt((4*A_nacl^3)/(27*B_nacl)))

[M_nacl, I_nacl] = max(Seq_nacl);

%%% Kappa value %%%
Kappa_nacl = (4 * A_nacl^3)/(27 * dp_nacl^3 * (log(Scrit_nacl))^2)




%% Plotting


semilogx(Dp,Seq_nh4so4, 'LineWidth',2.0)
hold on
semilogx(Dp,Seq_nacl,'LineWidth',2.0)
plot(Dp(I_nh4so4),M_nh4so4, '*', 'MarkerSize', 20)
plot(Dp(I_nacl),M_nacl, '*', 'MarkerSize', 20)
hold off
title('K{\"o}hler Curves', 'Interpreter','latex')
xlabel('Wet diameter, [m]', 'Interpreter','latex')
ylabel('$S_{eq}$', 'Interpreter','latex')
xlim([dp 1e-5])
ylim([0.98 1.011])
legend('Ammonium Sulfiade, $(NH_4)_2SO_4$', 'Sodium Chloride, $NaCl$', '$(NH_4)_2SO_4$ Critical Saturation Ratio', '$NaCl$ Critical Saturation Ratio', 'Interpreter', 'latex')





