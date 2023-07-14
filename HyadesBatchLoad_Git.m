%Analyse multiple Hyades .cdf, where the file names go from 'File1'
%to 'FileN'. Should have an accompanying Matlab file 'variables'. Saves the
%result in the original directory, and in a central repository (the paths
%for this will need to be changed at the bottom of the file)

clear all
%Ask user to select files to analyze
[FileName,PathName]=uigetfile('*.cdf','Select the INPUT DATA FILE(s)','MultiSelect','on');
[FileName2,PathName2]=uigetfile([PathName, '*.mat'],'Select the Variables File','MultiSelect','on');
%Goes through the filenames and reads the number, to record the valid files
%in the dataset.
for i=1:length(FileName)
    str = regexprep(char(FileName(i)),'.cdf','');
    str = regexprep(str,'File','');
    ValidFiles(i) = str2num(str);
end
  
%Generate empty arrays to store data for valid files in
IFARAll = [];
NeutronsAll =  [];
ImplosionVelocityAll = [];
OlsonCRAll = [];
GainAll = [];
EnergyAll = [];
TimeAll = [];
BurnUpAll = [];
MaxHSTempAll = [];
MaxHSRhoRAll = [];
MaxShellTempAll = [];
MaxShellRhoRAll = [];
RequiredLaserEnergyAll = [];
ErrorFiles = [];
HydroEfficiencyAll = [];
ShellKineticEnergyAll = [];
TotalEfficiencyAll = [];
VapourRadiusAll = [];
IceRadiusAll = [];
CHRadiusAll = [];
VapourPressureAll = [];
IceDensityAll = [];
CHDensityAll = [];


%Loop over all valid files, run the analysis script and save the key
%variables. Some files are missing, and some will error - these are not
%included in array. As such, the array position no longer corresponds to
%file number - ValidFiles contains the file number.
for FileIndex=ValidFiles
     
    try
file=[char(PathName) 'File' num2str(FileIndex) '.cdf']

Type=1;
Plots=0;

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


%A 0.8 factor is applied to input power to account for cross beam energy
%transfer (i.e., a power E/0.8 is required in real life to achieve E in
%Hyades. Laser Energy and Laser Power corresponds to real energy, while
%Simulation Power describes Hyades power.
LaserEnergy  = ncread(file,'Elasin')*(10^-7)/0.8;
LaserPower=[0; diff(LaserEnergy)./diff(Time)];
SimulationPower = LaserPower*0.8;

%Laser statistics calculated
TotalLaserEnergy = LaserEnergy(end);
TotalLaserPower = LaserPower;

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
Analysis_Git

if isempty(HydroEfficiencyValue)
    HydroEfficiencyValue = 0;
    ShellKineticEnergyValue = 0;
    TotalEfficiencyValue = 0;
    MaxOlsonCR = 0;
end

if isempty(MaxOlsonCR)
      MaxOlsonCR = 0;
end
IFARAll = [IFARAll, IFARCraxton(end)];
NeutronsAll = [NeutronsAll, max(Neutrons)];
ImplosionVelocityAll = [ImplosionVelocityAll, abs(min(ImplosionVelocity))];
OlsonCRAll = [OlsonCRAll, MaxOlsonCR];
GainAll = [GainAll, Gain];
EnergyAll = [EnergyAll, RequiredLaserEnergy];
TimeAll = [TimeAll, Time(end)];
BurnUpAll = [BurnUpAll, BurnUpFraction];
MaxHSTempAll = [MaxHSTempAll, PeakHSTemp];
MaxHSRhoRAll = [MaxHSRhoRAll, PeakHSRhoR];
MaxShellTempAll = [MaxShellTempAll, PeakShellTemp];
MaxShellRhoRAll = [MaxShellRhoRAll, PeakShellRhoR];
RequiredLaserEnergyAll = [RequiredLaserEnergyAll, RequiredLaserEnergy];
ShellKineticEnergyAll = [ShellKineticEnergyAll, ShellKineticEnergyValue];
HydroEfficiencyAll = [HydroEfficiencyAll, HydroEfficiencyValue];
TotalEfficiencyAll = [TotalEfficiencyAll, TotalEfficiencyValue];
VapourRadiusAll = [VapourRadiusAll, IceBoundary(1)];
IceRadiusAll = [IceRadiusAll, CHBoundary(1)];
CHRadiusAll = [CHRadiusAll, InitialRadius];
VapourPressureAll = [VapourPressureAll, Density(1,1)];
IceDensityAll = [IceDensityAll, Density(IceBoundaryIndex,1)];
CHDensityAll = [CHDensityAll, Density(CHBoundaryIndex,1)];

    catch 
        ErrorFiles = [ErrorFiles, FileIndex];
        ValidFiles(ValidFiles==FileIndex)=[];
    end


end

%Load Variables 
load([PathName2 FileName2])
PulseAll = PulseOnTime(:, ValidFiles);
PowerAll = PulsePower(:, ValidFiles);

if exist('HeatingTimings','var')==1
HeatingTimingsAll = HeatingTimings(ValidFiles, :);
HeatingAmpAll = Heating(ValidFiles, :);
end

if exist('HeatingEnergy','var')==1
HeatingEnergyAll = HeatingEnergy(ValidFiles);
end


if exist('SecondLaserOnTime','var')==1
    if length(SecondLaserOnTime)>1
SecondLaserOnTimeAll = SecondLaserOnTime(ValidFiles);
SecondLaserOffTimeAll = SecondLaserOffTime(ValidFiles);
SecondLaserPowerAll = SecondLaserPower(ValidFiles);
    else
SecondLaserOnTimeAll = repmat(SecondLaserOnTime,1,length(ValidFiles));
SecondLaserOffTimeAll = repmat(SecondLaserOffTime,1,length(ValidFiles));
SecondLaserPowerAll = repmat(SecondLaserPower,1,length(ValidFiles));
    end
end


%Delete unnecessary data from the import
clearvars -except ValidFiles IFARAll NeutronsAll ImplosionVelocityAll ...
    OlsonCRAll GainAll EnergyAll TimeAll BurnUpAll MaxHSTempAll ...
    MaxHSRhoRAll MaxShellTempAll MaxShellRhoRAll RequiredLaserEnergyAll ...
    RequiredLaserEnergyAll ErrorFiles FileName PathName HydroEfficiencyAll ...
    TotalEfficiencyAll ShellKineticAll VapourRadiusAll IceRadiusAll CHRadiusAll ... 
    VapourPressureAll IceDensityAll CHDensityAll PowerAll PulseAll ...
    ShellKineticEnergyAll Pulse1Time Pulse1TimeValues Pulse2Time ... 
    Pulse2TimeValues Pulse3Time Pulse3TimeValues Pulse4Time Pulse4TimeValues...
    RadiusVapour RadiusVapourValues RadiusIce RadiusIceValues RadiusCH RadiusCHValues...
    HeatingAmpAll HeatingTimingsAll HeatingEnergyAll PathName2 FileName2 ...
    SecondLaserOnTimeAll SecondLaserOffTimeAll SecondLaserPowerAll SecondLaserWavelength


%Save the data to a workspace in the directory, and back up on the H drive.
prompt = 'Save Data? y/n : ';
str = input(prompt,'s');
prompt = '4 pulses? y/n : ';
str2 = input(prompt,'s');
prompt = 'Auxilliary Heating? y/n : ';
str3 = input(prompt,'s');
prompt = 'Foam Only? y/n : ';
str4 = input(prompt,'s');
prompt = 'Second Laser? y/n: ';
str5 = input(prompt, 's');
prompt = 'Not yet backed up? y/n: ';
str6 = input(prompt, 's');

if str=='y' && str2=='n' && str3=='n' && str5=='n'
save([PathName, 'MatlabWorkspaceNew'])
BackUpName = strrep(PathName,'\','_');
BackUpName = strrep(BackUpName,':','');
save(['\\alfs1.physics.ox.ac.uk\al\paddock\Hyades Data\', BackUpName, '.mat'])
end


if str=='y' && str2=='y' && str3=='n' && str5=='n'
save([PathName, 'MatlabWorkspaceNew'])
BackUpName = strrep(PathName,'\','_');
BackUpName = strrep(BackUpName,':','');
save(['\\alfs1.physics.ox.ac.uk\al\paddock\Hyades Data\4 Pulses\', BackUpName, '.mat'])
end

if str3=='y' && str5=='n'
save([PathName, 'MatlabWorkspaceNew'])
BackUpName = strrep(PathName,'\','_');
BackUpName = strrep(BackUpName,':','');
save(['\\alfs1.physics.ox.ac.uk\al\paddock\Hyades Data\Auxilliary Heating\', BackUpName, '.mat'])
end

if str=='y' && str4=='y' && str5=='n'
save([PathName, 'MatlabWorkspaceNew'])
BackUpName = strrep(PathName,'\','_');
BackUpName = strrep(BackUpName,':','');
save(['\\alfs1.physics.ox.ac.uk\al\paddock\Hyades Data\Foam Only\', BackUpName, '.mat'])
end

if str5=='y'
save([PathName, 'MatlabWorkspaceNew'])
BackUpName = strrep(PathName,'\','_');
BackUpName = strrep(BackUpName,':','');
save(['\\alfs1.physics.ox.ac.uk\al\paddock\Hyades Data\Second Laser\', BackUpName, '.mat'])
end

if str6=='y'
save(['\\alfs1.physics.ox.ac.uk\al\paddock\Hyades Data\Not Yet Backed Up\', BackUpName, '.mat'])
end
