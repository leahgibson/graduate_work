% activation of populations of particles

%%% Each of these three parcels are involved in the development of a cloud
%%% that reaches a maximum supersaturation of 0.2%



clc;
clear all;
close all;

%% general parameters

surface_tension_water = 0.073; % J/m^2
Mw = 0.018; % molecular weight of wtaer, kg/mol
density_water = 1000; % kg/m^3
T = 298; % temperature, K
R = 8.314; % universal gas constant, J/molK

%% create bins

%%%% Dp space %%%%
start_val = 1; % units: nm
num_bins = 40;

% find upper value of each bin
bin_upper = zeros([1 num_bins]);
for i = 1:num_bins
    bin_upper(i) = start_val*(2^(1/3))^i;
end

% 41 edges for diameter of 1nm to 10.3 um
bin_edges_Dp = [start_val, bin_upper];

% mean diameter
mid_bin_diameter = zeros([1 num_bins]);
for i = 1:num_bins
    mid_bin_diameter(i) = sqrt(bin_edges_Dp(i)*bin_edges_Dp(i+1));
end

%% parcels

% parcle 1: sulfuric acid
sulfuric_N = 1000; % cm^-3
sulfuric_Dpg = 80; % nm
sulfuric_sigma = 1.6; 
sulfuric_Kappa = 0.9;

% parcel 2: organics
organics_N = 1100; % cm^-3
organics_Dpg = 110; % nm
organics_sigma = 1.8;
organics_Kappa = 0.1;

% parcel 3: mixture, 50% ammonium sulfate, 50% organics
mixture_N = 950; %cm^-3
mixture_Dpg = 90; %nm
mixture_sigma = 2.1;
mixture_Kappa = 0.5*organics_Kappa + 0.5*sulfuric_Kappa; 




%% plotting lognormal distributions in dN/dlogDp

n = length(mid_bin_diameter);


sulfuric_distribution = zeros([1 n]);
organics_distrbution = zeros([1 n]);
mixture_distribution = zeros([1 n]);
total_distribution = zeros([1 n]);

for i = 1:n
    sulfuric_distribution(i) = (log(10)*sulfuric_N/((2*pi)^(0.5) * log(sulfuric_sigma))) * ...
        exp(-(log(mid_bin_diameter(i)) - log(sulfuric_Dpg))^2 / (2*(log(sulfuric_sigma))^2));

    organics_distrbution(i) = (log(10)*organics_N/((2*pi)^(0.5) * log(organics_sigma))) * ...
        exp(-(log(mid_bin_diameter(i)) - log(organics_Dpg))^2 / (2*(log(organics_sigma))^2));

    mixture_distribution(i) = (log(10)*mixture_N/((2*pi)^(0.5) * log(mixture_sigma))) * ...
        exp(-(log(mid_bin_diameter(i)) - log(mixture_Dpg))^2 / (2*(log(mixture_sigma))^2));

    total_distribution(i) = sulfuric_distribution(i) + organics_distrbution(i) + mixture_distribution(i);
end


figure(1)
semilogx(mid_bin_diameter, sulfuric_distribution)
hold on
semilogx(mid_bin_diameter, organics_distrbution)
semilogx(mid_bin_diameter, mixture_distribution)
hold off
title('Distrubition of each Parcel', 'Interpreter', 'latex')
xlabel('Dry Diameter [nm]', 'Interpreter', 'latex')
ylabel('$\frac{dN}{d \log D_p}$ [cm$^{-3}$]', 'Interpreter', 'latex')
legend('Sulfuric acid', 'Organics', 'Ammonium sulfate and organics mixture')

figure(2)
semilogx(mid_bin_diameter,total_distribution)
title('Air Parcels Distribution', 'Interpreter', 'latex')
xlabel('Dry Diameter [nm]', 'Interpreter', 'latex')
ylabel('$\frac{dN}{d \log D_p}$ [cm$^{-3}$]', 'Interpreter', 'latex')

%% calculations for Scrit

A = (4 * Mw * surface_tension_water)/(R * T * density_water); % m
A = A*1e9; % nm

sulfuric_scrit = zeros([1 n]);
organics_scrit = zeros([1 n]);
mixture_scrit = zeros([1 n]);

for i = 1:n
    sulfuric_scrit(i) = exp(sqrt((4*A^3)/(27*mid_bin_diameter(i)^3*sulfuric_Kappa)));
    organics_scrit(i) = exp(sqrt((4*A^3)/(27*mid_bin_diameter(i)^3*organics_Kappa)));
    mixture_scrit(i) = exp(sqrt((4*A^3)/(27*mid_bin_diameter(i)^3*mixture_Kappa)));
end

figure(3)
plot(sulfuric_scrit, sulfuric_distribution)
hold on
plot(organics_scrit, organics_distrbution)
plot(mixture_scrit, mixture_distribution)
xlim([1 1.1])
hold off
xlabel('$S_{crit}$', 'Interpreter', 'latex')
ylabel('$\frac{dN}{d \log D_p}$ [cm$^{-3}$]', 'Interpreter', 'latex')
legend('Sulfuric acid', 'Organics', 'Ammonium sulfate and organics mixture')

figure(4)
plot(sulfuric_scrit, sulfuric_distribution)
hold on
plot(organics_scrit, organics_distrbution)
plot(mixture_scrit, mixture_distribution)
plot(1.002,0, 'k*', 'MarkerSize', 20)
xlim([1 1.01])
hold off
xlabel('$S_{crit}$', 'Interpreter', 'latex')
ylabel('$\frac{dN}{d \log D_p}$ [cm$^{-3}$]', 'Interpreter', 'latex')
legend('Sulfuric acid', 'Organics', 'Ammonium sulfate and organics mixture')

%% activation stats

max_supersat = 1.002;

sulfuric_concentration = 0;
organics_concentration = 0;
mixture_concentration = 0;

sulfuric_concentration_total = 0;
organics_concentration_total = 0;
mixture_concentration_total = 0;

% use midpoint riemann sums 
for i = 1:n
    if sulfuric_scrit(i) <= max_supersat
        sulfuric_concentration = sulfuric_concentration + (sulfuric_distribution(i)*(bin_edges_Dp(i+1)-bin_edges_Dp(i)));
    end

    if organics_scrit(i) <= max_supersat
        organics_concentration = organics_concentration + (organics_distrbution(i)*(bin_edges_Dp(i+1)-bin_edges_Dp(i)));
    end

    if mixture_scrit(i) <= max_supersat
        mixture_concentration = mixture_concentration + (mixture_distribution(i)*(bin_edges_Dp(i+1)-bin_edges_Dp(i)));
    end

    sulfuric_concentration_total = sulfuric_concentration_total + (sulfuric_distribution(i)*(bin_edges_Dp(i+1)-bin_edges_Dp(i)));
    organics_concentration_total = organics_concentration_total + (organics_distrbution(i)*(bin_edges_Dp(i+1)-bin_edges_Dp(i)));
    mixture_concentration_total = mixture_concentration_total + (mixture_distribution(i)*(bin_edges_Dp(i+1)-bin_edges_Dp(i)));
end

sulfuric_concentration
organics_concentration
mixture_concentration

sulfuric_concentration_total
organics_concentration_total
mixture_concentration_total

sulfuric_fraction = sulfuric_concentration/sulfuric_concentration_total
organics_fraction = organics_concentration/organics_concentration_total
mixture_fraction = mixture_concentration/mixture_concentration_total

%% diameter of particle at activation (i.e. when Scrit = 1.002)

sulfuric_diff = zeros([1 n]);
organics_diff = zeros([1 n]);
mixture_diff = zeros([1 n]);

for i = 1:n
    sulfuric_diff(i) = abs(sulfuric_scrit(i) - max_supersat);
    organics_diff(i) = abs(organics_scrit(i) - max_supersat);
    mixture_diff(i) = abs(mixture_scrit(i) - max_supersat);
end

[sulfuric_M,sulfuric_I] = min(sulfuric_diff);
[organics_M, organics_I] = min(organics_diff);
[mixture_m, mixture_I] = min(mixture_diff);

% finding clodest diameter to activation

sulfuric_diam = mid_bin_diameter(sulfuric_I)
organics_diam = mid_bin_diameter(organics_I)
mixture_diam = mid_bin_diameter(mixture_I)






