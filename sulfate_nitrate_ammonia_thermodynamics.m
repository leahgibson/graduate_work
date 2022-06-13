% Sulfate-Nitrate-Ammonia Thermodynamics

%%% This code plots the dynamics of the sulfate-nitrate-ammonia system. We
%%% assume the relative humidity is low enough that any ammonium nitrarte
%%% that forms is solid. 


clc;
close all;
clear all;

%% parameters and constants


R = 8.3145; % universal gas constant, J/mol K
molecular_weight_nh3 = 17/1000; % kg/mol
molecular_weight_nh4 = 18/1000; % kg/mol
molecular_weight_so4 = 96/1000; % kg/mol
molecular_weight_nh42so4 = 132/1000; % kg/mol
molecular_weight_h2so4 = 98/1000; % kg/mol
molecular_weight_hno3 = 63/1000; % kg/mol
molecular_weight_nh4no3 = 80/1000; % kg/mol
pressure = 100000; % pascals

T = 280; % temperature in Kelvin
total_hno3_ppbv = 15; % ppbv
total_nh3_ppbv = 20; % ppbv



%% set up sulfate

% sulfate concentrations
n = 1000;
max_sulfate = 20; % ppbv
max_sulfate_molar_conc = (max_sulfate*1e-9*pressure)/(R*T); % mol/m^3
total_sulfate = linspace(0,max_sulfate_molar_conc,n); % mol/m^3




%% compute free ammonia


% convert ppbv to molar concentration
total_ammonia = (total_nh3_ppbv*1e-9*pressure)/(R*T); % mol/m^3

free_ammonia = zeros([1 n]);
ammonium_sulfate_particles_mass = zeros([1 n]);
free_ammonia_ppbv = zeros([1 n]);
for i = 1:n
    free_ammonia(i) = total_ammonia-2*total_sulfate(i); % mol/m^3
    free_ammonia_ppbv(i) = (free_ammonia(i)*R*T)/(pressure*1e-9);
    % amount of annomium sulfate formed, ug/m^3
    if free_ammonia(i) >= 0 % only formed when there is enough ammonia
        ammonium_sulfate_particles_mass(i) = total_sulfate(i)*molecular_weight_nh42so4*1e9; 
        I = i; % will save index of largest ammount of ammonium sulfate formed before annomia ran out
    else
        remaining_sulfate = total_sulfate(i)-total_sulfate(I); % left over sulfate that will be ammonium sulfate
        ammonium_sulfate_particles_mass(i) = remaining_sulfate*molecular_weight_h2so4*1e9 + ammonium_sulfate_particles_mass(I);
    end
end

%% gas ratio and equilibrium constant

% convert ppbv nitrate to molar concentration
total_nitrate = (total_hno3_ppbv*1e-9*pressure)/(R*T); % mol/m^3

% gas ratio
gas_ratio = zeros([1 n]);

for i = 1:n
    gas_ratio(i) = free_ammonia(i)/total_nitrate;
end

% equilibrium constant, Kp in ppbv^2
Kp = exp(84.6-(24220/T)-6.1*log(T/289));

%% formation of ammonium nitrate

ammonium_nitrate_particles_mass = zeros([1 n]);
total_particle_mass = zeros([1 n]);
ammonium_nitrate_concentration = zeros([1 n]);

for i = 1:n

    if gas_ratio(i) < 0
        disp('no ammonium nitrate formed for sulfate concentration of ' + string(total_sulfate(i)))
        total_particle_mass(i) = ammonium_sulfate_particles_mass(i);

    elseif gas_ratio(i) >= 0 && gas_ratio(i) <= 1
        disp('ammonia limited for sulfate concentration of ' + string(total_sulfate(i)))
        if free_ammonia_ppbv(i)*total_hno3_ppbv >= Kp
            r = roots([1 -(total_hno3_ppbv + free_ammonia_ppbv(i)) (free_ammonia_ppbv(i)*total_hno3_ppbv - Kp)]);
            solution = r(r < total_hno3_ppbv & r < free_ammonia_ppbv(i) & r > 0);
            if size(solution) > 1
                disp('error')
            end
            ammonium_nitrate_concentration(i) = solution; % ppbv
            ammonium_nitrate_particles_mass(i) = (ammonium_nitrate_concentration(i)*pressure*molecular_weight_nh4no3)/(R*T); % ug/m^3
            total_particle_mass(i) = ammonium_sulfate_particles_mass(i) + ammonium_nitrate_particles_mass(i);
        else
            disp('no ammonium nitrate formed for sulfate concentration of ' + string(total_sulfate(i)))
            total_particle_mass(i) = ammonium_sulfate_particles_mass(i);
        end

    elseif gas_ratio(i) > 1
        disp('nitric acid limited for sulfate concentration of ' + string(total_sulfate(i)))
        if free_ammonia_ppbv(i)*total_hno3_ppbv >= Kp
            r = roots([1 -(total_hno3_ppbv + free_ammonia_ppbv(i)) (free_ammonia_ppbv(i)*total_hno3_ppbv - Kp)]);
            solution = r(r < total_hno3_ppbv & r < free_ammonia_ppbv(i) & r > 0);
            if size(solution) > 1
                disp('error')
            end
            ammonium_nitrate_concentration(i) = solution; % ppbv
            ammonium_nitrate_particles_mass(i) = (ammonium_nitrate_concentration(i)*pressure*molecular_weight_nh4no3)/(R*T); % ug/m^3
            total_particle_mass(i) = ammonium_sulfate_particles_mass(i) + ammonium_nitrate_particles_mass(i);
        else
            disp('no ammonium nitrate formed for sulfate concentration of ' + string(total_sulfate(i)))
            total_particle_mass(i) = ammonium_sulfate_particles_mass(i);
        end 
    end
end


%% no nitric acid senario

total_hno3_ppbv = 15; % ppbv



%% figures

figure(1)
plot(total_sulfate*molecular_weight_so4*1e9, total_particle_mass, 'Linewidth', 2)
hold on
plot(total_sulfate*molecular_weight_so4*1e9, ammonium_sulfate_particles_mass, 'Linewidth', 2)
hold off
legend('Nitric Acid = 15 ppbv', 'Nitric acid = 0 ppbv')
title('Sulfate-Nitrate-Ammonia Thermodynamics')
xlabel('$SO_4^{2^-} [\mu g/m^{-3}]$', 'Interpreter','latex')
ylabel('Total particle mass, $[\mu g/m^{-3}]$', 'Interpreter','latex')











