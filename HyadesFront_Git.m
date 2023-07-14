%This script takes a .cdf Hyades Output and imports key variables. It also
%creates a few other key variables from this data. Once the imported data
%is in the correct form, it is sent to the Analysis script. NOTE THAT
%ADIABAT IS INCORRECT AND CANNOT BE TRUSTED.

%Insert filepath. ncdisp will give info about the variables contained
%within the file that can be taken


clear all

file = 'C:\Users\paddock\Documents\Hyades\SavedFiles\Cutoffs\0.25Size\0.25Size4PulseCutoff.cdf'

% ncdisp(file) %This command allows you to see all saved variables (uncomment to use)

%I have 4 modes for how many plots will open - 0 is none, 1 is just an
%output summary page, 1 also includes some important plots to check that
%the analysis script has worked possibly, 2 also includes a density plot,
%and 3 are lots of plots that I've used for different purposes.
Plots=3;

%Import the key variables, and get in correct format
Radius  = ncread(file,'R');
Time  = ncread(file,'DumpTimes');
Radiuscm  = ncread(file,'Rcm');
Density  = ncread(file,'Rho'); %g/cm^3
ElecTemp  = ncread(file,'Te')*1000; %Converted to eV from keV
IonTemp  = ncread(file,'Ti')*1000; %Converted to eV from keV
Neutrons = sum((0.8/(14.1 * 1.60218e-6))*ncread(file, 'Bpeprd')); %x 0.8 to get fraction of energy in neutrons, then divide by neutron energy in MeV
Neutrons(isnan(Neutrons)) = 0;
NeutronRate = sum((0.8/(14.1 * 1.60218e-6))*ncread(file, 'Bpeprdr'));
OldRhoR = Density.*Radius(1:end-1, :); 
Volume =  ncread(file,'Vol'); %cm^3
Mass = Volume(2:end,:).*Density;
Pressure = ncread(file,'Pres').*0.0000001; %Converted to J/cm^3 from dyn/cm^2
PressureGPa = Pressure.*10^-3;
TNOutput = ncread(file, 'Bpeprdr'); %TN burn rate
TNInput = ncread(file, 'Bpedep').*10^-7; %Deposited TN energy per zone, converted from erg to J
TNInputRate = [zeros(size(TNInput,1),1),diff(TNInput,1,2)]; %TN deposition rate
Velocity = ncread(file, 'U')./100000; %Converted to km/s from cm/s
IonThermalEnergy  = ncread(file,'Eion').*10^-7; %Ion Thermal Energy, converted from erg to J.
ElectronThermalEnergy  = ncread(file,'Eelc').*10^-7; %Electron Thermal Energy, converted from erg to J.
%KineticEnergy  = ncread(file,'Ekint').*10^-7;
DepositedLaserPowerZones = ncread(file,'Deplas').*10^-7; %Deposited Laser Energy per zone, converted from erg to J.
DepositedLaserPower = sum(DepositedLaserPowerZones).'; %Total deposited laser power (over all zones)
DepositedEnergy = trapz(Time, DepositedLaserPower); %Integrate for deposited energy
try
    zbar = ncread(file,'Zbar');
end

%A 0.8 factor is applied to input power to account for cross beam energy
%transfer (i.e., a power E/0.8 is required in real life to achieve E in
%Hyades). Laser Energy and Laser Power corresponds to real energy, while
%Simulation Power describes Hyades power.
LaserEnergy  = ncread(file,'Elasin')*(10^-7)/0.8;
LaserPower=[0; diff(LaserEnergy)./diff(Time)];
SimulationPower = LaserPower*0.8;

%Laser statistics calculated
TotalLaserEnergy = LaserEnergy(end); %Laser energy used over whole simulation
TotalLaserPower = LaserPower;



%For the output summary, I have included an indication of the different
%pulse powers and when power changes. The following two chunks of code
%estimate this:

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




%Run Analysis script now that data is in correct format (the 'type' just tells
%the analysis script that this is a Hyades output, as I've also used other
%codes in the past - these use a seperate front end, but run the same analysis code)
Type=1;
Analysis_Git

