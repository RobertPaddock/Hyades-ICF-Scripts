%This script takes a .exo Helios Output and imports key variables. It
%operates similary to HyadesFront - the Helios data is pre-processed, so that
%the same analysis script can be used.

%Insert filepath. ncdisp will give info about the variables contained
%within the file that can be taken
file = '\\aldaq1.physics.ox.ac.uk\Archer\Robert\HeliosSims\080223\IgorDT\Mixd_Igor_1pt35.exo'




%ncdisp(file); 

Plots = 1;

%Below are some key files
%2 shell change second pulse time
%ReducedPower2Amp_4xPressure
%4xSize_IncreasedAR8_2
%2ndLaserEarly
%File13IncreasedResolution.exo
%File43IncreasedResolution - Compare with Hyades

%Import data from file. Ice and CH Boundary depend on the Vapour having 40
%zones and the Ice having 80.
Time  = ncread(file,'time_whole');
Radius  = ncread(file,'zone_boundaries');
Density  = ncread(file,'mass_density');
IonTemp  = ncread(file,'ion_temperature');
OldRhoR = Density.*Radius(1:end-1, :); 
try
Neutrons=ncread(file, 'TimeIntFusionProd_n_1406');
Neutrons=Neutrons.';
catch
    Neutrons = zeros(1, length(Time));
end
NeutronRate = diff(Neutrons);

Mass=ncread(file,'zone_mass');
Pressure=ncread(file, 'ion_pressure'); %in J/ cm^3
PressureGPa = Pressure.*10^-3;
Velocity=ncread(file, 'fluid_velocity')./100000; %Converted to km/s from cm/s
Volume = Mass./Density;


%Helios needs a factor 0.6 correction to energy (i.e. E/0.6 is required to
%deliver E of energy in Helios). Laser Energy/Laser Power is real,
%experimental energy and power, while Simulation Power is the value used
%by Helios.
LaserEnergy=ncread(file,'LaserEnDeliveredTimeInt')/0.8;
%SimulationPower = ncread(file,'LaserPwrDeliveredForBeam');
%LaserPower = SimulationPower.'/0.6;
LaserPower=[0; diff(LaserEnergy)./diff(Time)];
SimulationPower = LaserPower*0.8;

%Laser statistics calculated
TotalLaserEnergy = LaserEnergy(end);
TotalLaserPower = LaserPower;

DepositedLaserPowerZones = ncread(file,'LaserPwrSrc'); %Laser energy deposition rates [J per (g.sec)].
DepositedLaserPower = sum(DepositedLaserPowerZones.*Mass);
DepositedEnergy1 = trapz(Time, DepositedLaserPower);
DepositedEnergy = ncread(file,'LaserEnTimeIntg'); %Time-integrated laser energy deposition [J per (cm**X.zone)] .
DepositedEnergy = ncread(file,'LaserEnTimeIntg').*Volume; %Time-integrated laser energy deposition [J per (cm**X.zone)] .
DepositedEnergy = sum(DepositedEnergy(:, end));

%Identify where power changes. Use this to find the first index for each
%plateau, and then use this to find the pulse values.
LaserPowerDiff=(([0; diff(LaserPower)]./LaserPower)<0.0001).';
PulseIndices = findstr([0 LaserPowerDiff], [0 1])-1;
PulsePowers = LaserPower(PulseIndices);

%Using these indices and powers, calulate y=mx+c for each (by referring to
%the two points before the plateau). Then, use this equation to find the
%pulse turn on/off times. The subtraction at the end is a fudge factor.
%THis is only correct to the accuacy of the timesteps.
m = (LaserPower(PulseIndices-1) - LaserPower(PulseIndices-2))./(Time(PulseIndices-1) - Time(PulseIndices-2));
c = LaserPower(PulseIndices-1) - m.*Time(PulseIndices-1);
PulseOnTime = ((LaserPower([1 PulseIndices(1:end-1)])-c)./m) - 0.0005E-7;
PulseOnAndRiseTime = (LaserPower([PulseIndices(1:end)])-c(1:end))./m(1:end)- 0.0005E-7;
%num2str((abs(PulsePowers)./10^12).', '%.2f TW\n');

%Run Analysis script now that data is in correct format
Type = 0;
AnalysisCondensed

%Stagnation plot to show where CR was recorded
figure;
xlim([xlowerlim xupperlim]);
title('Convergence Ratio')
xlabel('Time (ns)');
ylabel('Radius (cm)');
hold on
plot(Time,OlsonCR, 'linestyle', '--','color', 'k');
xline(BangTime, 'r');
xline(CraxtonStagnation, 'b');
plot(Time(MaxOlsonCRIndex),OlsonCR(MaxOlsonCRIndex) ,'r*')
legend({'Convergence Ratio', 'Bang Time', 'Stagnation'}, 'Location', 'NorthWest');
hold off

max(max(Pressure))/10^8
max(max(Density))