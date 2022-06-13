% Organic Aerosol partitoning and aging of wood smoke aerosol

%%% This code tracks the total particle mass of smoke emitted from a
%%% wildfire (could also apply to stove emissions with different
%%% concentrations and dilution timescales)

clc;
clear all;
close all;

%% initial parameters

c_star = [.01, 0.1, 1, 10, 100, 1000, 10000]; % pure species equilibrium vapor concentration, ug/m^3
f_emissions = [0.2, 0, 0.1, 0.1, 0.2, 0.2, 0.3]; % total organic emissions
n = length(c_star);
total_organics = 5000; % ug/m^3

%% concentration of gas+particles in each bin

c = zeros([1 n]);
for i = 1:n
    c(i) = f_emissions(i)*total_organics; % ug/m^3
end

%% computing c_oa and f(i)

%C_OA := total mass concentration of particle phase organics
%f(i) := fraction in aerosol phase

aerosol_fraction = zeros([1 n]); % f(i)'s
particle_mass_concentration_guess = 5; % ug/m^3

for i = 1:n
    aerosol_fraction(i) = 1/(1 + (c_star(i)/particle_mass_concentration_guess));
end

particle_mass_concentration = sum(aerosol_fraction(1:n).*c(1:n));

while abs(particle_mass_concentration-particle_mass_concentration_guess) > 0.001 
    particle_mass_concentration_guess = particle_mass_concentration;
    for i = 1:n
        aerosol_fraction(i) = 1/(1 + (c_star(i)/particle_mass_concentration_guess));
    end
    particle_mass_concentration = sum(aerosol_fraction(1:n).*c(1:n));
end

% final update
for i = 1:n
    aerosol_fraction(i) = (1 + (c_star(i)/particle_mass_concentration))^(-1);
end

%% compute amount of bin that is aerosol

aerosol_amount = zeros([1 n]);

for i = 1:n
    aerosol_amount(i) = c(i)*aerosol_fraction(i);
end

%% figures

figure(1)
bar(log10(c_star),c)
hold on
bar(log10(c_star),aerosol_amount, 'FaceColor',[0 0.7 0.7])
set(gca,'Xtick',-2:4); %// adjust manually; values in log scale
set(gca,'Xticklabel',10.^get(gca,'Xtick')); %// use labels with linear values
xlabel('C$^*$ [$\mu$g m$^{-3}$]', 'Interpreter', 'latex')
ylabel('Mass concentration [$\mu$g m$^{-3}$]', 'Interpreter','latex')
legend('Gas', 'Aerosol Particle', 'Location','northwest')
hold off

%% dilution

% calculate every half hour
steps = 49;
time = linspace(0, 48, steps);

% total concentration of gas+particles in each bin for each time step
c_new = zeros([steps n]);
for t = 1:steps
    for i = 1:n
        c_new(t,i) = c(i)*exp(-time(t)/12);
    end
end

% total concentration of gas+particles for each time step
total_concentration = zeros([1 steps]);
for i = 1:steps
    total_concentration(i) = sum(c_new(i,(1:n)));
end

% calculate new C_OA and f(i)'s for each bin and time
particle_mass_concentrations_time = zeros([ 1 steps]);
aerosol_fractions_time = zeros([steps n]);

for t = 1:steps

    particle_mass_concentration_guess = 5; % ug/m^3

    for i = 1:n
        aerosol_fractions_time(t,i) = 1/(1 + (c_star(i)/particle_mass_concentration_guess));
    end

    particle_mass_concentrations_time(t) = sum(aerosol_fractions_time(t,(1:n)).*c_new(t,(1:n)));

    while abs(particle_mass_concentrations_time(t)-particle_mass_concentration_guess) > 0.001
        particle_mass_concentration_guess = particle_mass_concentrations_time(t);
        for i = 1:n
            aerosol_fractions_time(t,i) = 1/(1 + (c_star(i)/particle_mass_concentration_guess));
        end
        particle_mass_concentrations_time(t) = sum(aerosol_fractions_time(t,(1:n)).*c_new(t,(1:n)));
    end

    % final update
    for i = 1:n
        aerosol_fractions_time(t,i) = 1/(1 + (c_star(i)/particle_mass_concentration_guess));
    end

end

% fraction of organics in particle phase, C_OA/C
organics_fraction = zeros([1 n]);
for i = 1:steps
    organics_fraction(i) = particle_mass_concentrations_time(i)/total_concentration(i);
end



%% figures 

figure(2)
plot(time,particle_mass_concentrations_time, 'LineWidth', 2)
xlabel('Time, [hr]', 'Interpreter','latex')
ylabel('Particle-phase organic mass','Interpreter','latex')


figure(3)
plot(time,organics_fraction, 'LineWidth', 2)
xlabel('Time, [hr]', 'Interpreter','latex')
ylabel('Fraction of total organics in particle phase, $C_{OA}/C$','Interpreter','latex')

figure(4)
bar(log10(c_star),c_new(steps,:))
hold on
bar(log10(c_star),c_new(steps,:).*aerosol_fractions_time(steps,:), 'FaceColor',[0 0.7 0.7])
set(gca,'Xtick',-2:4); %// adjust manually; values in log scale
set(gca,'Xticklabel',10.^get(gca,'Xtick')); %// use labels with linear values
xlabel('C$^*$ [$\mu$g m$^{-3}$]', 'Interpreter', 'latex')
ylabel('Mass concentration [$\mu$g m$^{-3}$]', 'Interpreter','latex')
title('Final Distribution')
legend('Gas', 'Aerosol Particle', 'Location','northwest')
hold off






