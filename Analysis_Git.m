%This script analyzes the data imported via the front scripts (such as HyadesFront).
%It is called within the HyadesFront script.

%Identify Ice and CH layers based on density changes.
IceBoundaryIndex=min(find(Density(:,1)>0.211));
CHBoundaryIndex=min(find(Density(:,1)>1));
IceBoundary=Radius(IceBoundaryIndex,:);
CHBoundary=Radius(CHBoundaryIndex,:);

%Define conditions for HS and shell thresholds (Both thresholds should have values of 1
%inside shell and values of 0 outside i.e. Shell threshold will have a 1 if
%the zone is within the shell outer boundary, but a 0 otherwise). The criteria is outlined in
%the paper, but is based on identifying the max density 
HSThreshold=Density-Density(1,:) > (max(Density-Density(1,:))*0.135);
ShellThreshold=Density>max(Density)*0.135;

%Use HS threshold condition to define the HS radius
[~, HSIndex] =max(HSThreshold);
HSIndexReshaped=sub2ind(size(Radius),HSIndex,1:length(Time));
HSRadius=Radius(HSIndexReshaped);

%Use shell threshold condition to define the shell radius
ShellIndex = arrayfun(@(x)find(ShellThreshold(:,x),1,'last'),1:size(ShellThreshold,2), 'UniformOutput',false);
tf = cellfun('isempty',ShellIndex); % true for empty cells
ShellIndex(tf) = {1};
ShellIndex= [ShellIndex{:}];
ShellIndexReshaped=sub2ind(size(Radius),ShellIndex,1:length(Time));
ShellRadius=Radius(ShellIndexReshaped);

%Calculate Isentrope Parameter (adiabat)
DegeneratePressure = (2.17*10^12 * 1e-7) *(Density.^(5/3));
IsentropeParameter = Pressure./DegeneratePressure;

%Create array of 1 if in hotspot and 0 outside, and equivalent for shell
%(unlike shell threshold, this is 0 if in hotspot, 1 if in shell, and 0
%outside shell)
Volume = Mass./Density;
HS = zeros(size(Density,1), length(Time));
Shell = zeros(size(Density,1), length(Time));
for i=1:length(Time)
    HS(1:HSIndex(i), i) = 1;
    Shell(HSIndex(i): ShellIndex(i), i)=1;
end

%Measure mass of vapour zones
VapourZones = zeros(size(Density,1), length(Time));
VapourZones(1:IceBoundaryIndex, :) = 1;
VapourMass = VapourZones.*Mass;
TotalVapourMass = sum(VapourMass);

%Create array recording relevant values in HS, and take averages.
HSVolume = HS.*Volume;
HSMass = HS.*Mass;
HSTemp = HS.*IonTemp;
HSDensity = HS.*Density;
HSOldRhoR = HS.*OldRhoR;
HSAverageTemp = sum(HSMass.*HSTemp)./sum(HSMass);
HSAverageDensity = sum(HSMass.*HSDensity)./sum(HSMass);
HSAverageOldRhoR = sum(HSMass.*HSOldRhoR)./sum(HSMass); %An old definition of RhoR
HSVelocity = HS.*Velocity(2:end, :);
try
HSKineticEnergy = 0.5*sum((HSMass/1000).*((1000*HSVelocity).^2));
HSThermalEnergy = sum(HS.*(IonThermalEnergy+ElectronThermalEnergy));
end

%Create array recording relevant values in Shell, and take averages.
ShellVolume = Shell.*Volume;
ShellVelocity = Shell.*Velocity(2:end, :);
ShellIsentrope = Shell.*IsentropeParameter;
ShellMass = Shell.*Mass;
ShellTemp = Shell.*IonTemp;
ShellDensity = Shell.*Density;
ShellOldRhoR = Shell.*OldRhoR;
ShellAverageTemp = sum(ShellMass.*ShellTemp)./sum(ShellMass);
ShellAverageDensity = sum(ShellMass.*ShellDensity)./sum(ShellMass);
ShellAverageOldRhoR = sum(ShellMass.*ShellOldRhoR)./sum(ShellMass); %An old definition of RhoR
ShellKineticEnergy = 0.5*sum((ShellMass/1000).*((1000*ShellVelocity).^2));
try
ShellThermalEnergy = sum(Shell.*(IonThermalEnergy+ElectronThermalEnergy));
ShellMomentum = sum((ShellMass/1000).*((1000*ShellVelocity)));
end

try
%Calculate Total Energies of HS, Shell and Capsule
TotalHSEnergy = HSKineticEnergy + HSThermalEnergy;
TotalShellEnergy = ShellKineticEnergy + ShellThermalEnergy;
TotalEnergy = TotalHSEnergy+TotalShellEnergy;
TotalThermalEnergy = HSThermalEnergy+ShellThermalEnergy;
TotalKineticEnergy = HSKineticEnergy+ShellKineticEnergy;
end

%Calculate convergence ratio
InitialRadius = max(Radius(:,1));
OlsonCR = IceBoundary(1)./HSRadius; %We use Olson's definition
LindlCR = InitialRadius./HSRadius;

%Calculate Parametric Criteria
ParametricLimit = (max(TotalLaserPower))*(0.35^2)/(4*pi()*(Radius(end,1)^2));

%Define time period around time of convergence/stagnation
%If neutrons are switched on, use this. Otherwise, look for peak
%convergence and look around there
if max(Neutrons)>0
    MinStagnationIndex = find(Neutrons>max(Neutrons)/100000, 1,'first');
    MaxStagnationIndex = find(Neutrons>max(Neutrons)*0.9999, 1,'first');
    RequiredLaserEnergy = LaserEnergy(MaxStagnationIndex);
    RequiredDepositedEnergy = trapz(Time(1:MaxStagnationIndex), DepositedLaserPower(1:MaxStagnationIndex));
    yupperlim = min((ShellRadius(MinStagnationIndex:MaxStagnationIndex)))*1.5;
else 
    MinStagnationIndex = find(HSRadius<3*min(HSRadius), 1,'first');
    MaxStagnationIndex = max(find(HSRadius<2*min(HSRadius)));
    RequiredLaserEnergy = LaserEnergy(MaxStagnationIndex);
    RequiredDepositedEnergy = trapz(Time(1:MaxStagnationIndex), DepositedLaserPower(1:MaxStagnationIndex));
    yupperlim = min(HSRadius)*1.5;
end
xlowerlim = Time(MinStagnationIndex);
xupperlim =  Time(MaxStagnationIndex);

%Calculate neutron energy using required Laser Energy
NeutronEnergy = 14.06e6*1.602E-19*max(Neutrons);
Gain = NeutronEnergy/RequiredLaserEnergy;

%Calculate IFAR, defined as ablation front radius divided by shell
%thickness, at the time where ablation front radius is at 2/3 of the
%initial inner radius of the shell
IFARIndex=find(ShellRadius<2*IceBoundary(1)/3, 1,'first');
IFAR = ShellRadius./(ShellRadius - HSRadius);
IFARCraxton = IFAR(IFARIndex);

%Calculate Implosion Velocity
ImplosionVelocity = (sum(ShellVolume.*ShellVelocity)./sum(ShellVolume));
 
%Calculate min and average shell isentrope
CHZones = Density(:,1)>1;
NonCHZones = repmat((ones(size(Density,1), 1) - CHZones), 1,length(Time));
%AveragedIsentrope = sum(ShellVolume.*ShellIsentrope.*NonCHZones)./sum(ShellVolume.*NonCHZones);
AveragedIsentrope = sum(ShellMass.*ShellIsentrope.*NonCHZones)./sum(ShellMass.*NonCHZones);
ShellIsentrope(ShellIsentrope==0)=NaN;
MinIsentrope = min(ShellIsentrope);

%Calculate Hydroefficiency, total efficiency, and find max values for these
%and shell Ek
%Commented lines look for max before shell starts moving in opposite
%direction, but this can give wrong result if time resolution is too
%coarse. New lines instead look for max in from when file starts to the
%first time after maximum velocity when the kinetic energy drops under 0.7
%times what it was at max velocity
HydroEfficiency = (ShellKineticEnergy./RequiredDepositedEnergy)*100;
TotalEfficiency = (ShellKineticEnergy./RequiredLaserEnergy)*100;
% HydroEfficiency = (ShellKineticEnergy./DepositedEnergy)*100;
% TotalEfficiency = (ShellKineticEnergy./TotalLaserEnergy)*100;
MaxVelocityIndex = find(ImplosionVelocity==min(ImplosionVelocity), 1, 'first');
ShellKineticEnergyValue = max(ShellKineticEnergy(1:(MaxVelocityIndex+find(ShellKineticEnergy(MaxVelocityIndex:end)<0.7*ShellKineticEnergy(MaxVelocityIndex), 1, 'first'))));
HydroEfficiencyValue = max(HydroEfficiency(1:(MaxVelocityIndex+find(ShellKineticEnergy(MaxVelocityIndex:end)<0.7*ShellKineticEnergy(MaxVelocityIndex), 1, 'first'))));
TotalEfficiencyValue = max(TotalEfficiency(1:(MaxVelocityIndex+find(ShellKineticEnergy(MaxVelocityIndex:end)<0.7*ShellKineticEnergy(MaxVelocityIndex), 1, 'first'))));
%ShellKineticEnergyValue = max(ShellKineticEnergy(1:find(sum(ShellVelocity)>0, 1, 'first')));
%HydroEfficiencyValue = max(HydroEfficiency(1:find(sum(ShellVelocity)>0, 1, 'first')));
%TotalEfficiencyValue = max(TotalEfficiency(1:find(sum(ShellVelocity)>0, 1, 'first')));




%Calculate the HS, Shell, and overall capsule RhoR using the new definition
%of RhoR (integrating density over radius). Do this at each timestep.
for i=1:length(Time)
    try
HSRhoR(i) = trapz(Radius(1:HSIndex(i), i), Density(1:HSIndex(i), i));
    catch 
        HSRhoR(i)=0;
    end
ShellRhoR(i) = trapz(Radius(HSIndex(i):ShellIndex(i), i), Density(HSIndex(i):ShellIndex(i), i));
CapsuleRhoR(i) = trapz(Radius(1:ShellIndex(i), i), Density(1:ShellIndex(i), i));

%For above we integrate up to shell boundary. Below checks that this isn't
%including CH ablator - if it is, only go up to the CH boundary, so that
%only DT is included in RhoR calculation.
if ShellIndex(i)<CHBoundaryIndex
CapsuleRhoRDTOnly(i) = trapz(Radius(1:ShellIndex(i), i), Density(1:ShellIndex(i), i));
else
    CapsuleRhoRDTOnly(i) = trapz(Radius(1:CHBoundaryIndex, i), Density(1:CHBoundaryIndex, i));
end

end

%Calculate max CR. Idenitify any major jumps in CR by looking for rapid
%changes in gradient. Otherwise, run till the end of the sim. Find the max
%CR in the relevant range. CR can be quite fiddly, so need to check this
CRWindowMin = MaxVelocityIndex;
% CRWindowMax = MaxVelocityIndex + 1 + find(diff(smooth(OlsonCR(MaxVelocityIndex:end),5))<0, 1, 'first');
CRWindowMax = MaxVelocityIndex + find(abs(diff(diff(OlsonCR(MaxVelocityIndex:end))))>2.5, 1, 'first') -1;
if isempty(CRWindowMax)
    CRWindowMax = length(Time);
end
[MaxOlsonCR, MaxOlsonCRIndex] = max(OlsonCR(CRWindowMin:CRWindowMax));
[MaxLindlCR, MaxLindlCRIndex] = max(LindlCR(CRWindowMin:CRWindowMax));
MaxOlsonCRIndex = CRWindowMin+MaxOlsonCRIndex-1;
% figure
% yyaxis left
% plot(OlsonCR)
% ylim([0 20])
% yyaxis right
% plot(diff(diff(OlsonCR)))
% ylim([-200 200])
% ylim([-10 10])

%Simulate some data, and plot the Ion Temp vs OldRhoR against Cheng's analytic
%formula
SimT = 4000:10:10000;
SimOldRhoR = 4*5.514*(SimT/1000).^-2.5;
TimestepsAboveIgnition=sum(HSAverageOldRhoR>4*5.514*(HSAverageTemp/1000).^-2.5);

%Quick estimate of Burn Up Fraction;
MassVapour = sum(Mass(1:IceBoundaryIndex, 1));
MassDT = sum(Mass(1:CHBoundaryIndex, 1));
MolesDT = MassDT/5;
MoleculesDT = MolesDT*6.022*10^23;
BurnUpFraction = (max(Neutrons)/MoleculesDT)*100;

%Burn averaged quantities, from Lindl, Physics of Plasmas, 25, 122704. 
TimeLength = diff(Time);
BurnAvTemp = sum((TimeLength.').*sum(IonTemp(:,2:end).*TNOutput(:,2:end)))./sum((TimeLength.').*sum(TNOutput(:,2:end)));
BurnAvPressure = sum((TimeLength.').*sum(PressureGPa(:,2:end).*TNOutput(:,2:end)))./sum((TimeLength(:,2:end).').*sum(TNOutput(:,2:end)));

%Find Bang Time and Stagnation Time (according to Craxton)
NeutronsPerTime = [0, diff(Neutrons)];
[~, BangTimeIndex] = max(NeutronsPerTime);
BangTime = Time(BangTimeIndex);
[~, CraxtonStagnationIndex] = max(sum(Pressure));
CraxtonStagnation = Time(CraxtonStagnationIndex);
BurnWidth = Time(find(NeutronRate>0.5*max(NeutronRate), 1, 'last')) - Time(find(NeutronRate>0.5*max(NeutronRate), 1, 'first'));

%Find max average T and OldRhoR for HS and Shell;
PeakHSRhoR = max(HSRhoR(MinStagnationIndex:MaxStagnationIndex));
PeakHSTemp = max(HSAverageTemp(MinStagnationIndex:MaxStagnationIndex));
PeakShellRhoR = max(ShellRhoR(MinStagnationIndex:MaxStagnationIndex));
PeakShellTemp = max(ShellAverageTemp(MinStagnationIndex:MaxStagnationIndex));

%Calculate alpha deposition (technically in whole capsule, but assume this
%is just HS)
if Type == 1
OverallTNInput = sum(TNInput);
OverallTNEnergyDeposited = OverallTNInput(end);
BangTimeTNEnergyDeposited = OverallTNInput(BangTimeIndex);
end

try
HSTNEnergyDeposited = HS.*TNInputRate;
TotalHSTNEnergyDepositionRate = sum(HSTNEnergyDeposited);
TotalHSTNEnergyDeposited = sum(TotalHSTNEnergyDepositionRate);
BangTimeHSTNEnergyDeposited = sum(TotalHSTNEnergyDepositionRate(1:BangTimeIndex));
end

% CHBoundary=IceBoundary
% CHBoundaryIndex=IceBoundaryIndex


try
%Identify time after max implosion velocity where kinetic and thermal
%energies equal one another. This is an approximate value for the point at
%which the energy is flat before alpha deposition.
[~,FlatEnergyIndex] = max(TotalThermalEnergy(MaxVelocityIndex:end) > TotalKineticEnergy(MaxVelocityIndex:end));
FlatEnergyIndex = MaxVelocityIndex+FlatEnergyIndex;
FlatEnergyTime = Time(FlatEnergyIndex);


%Use the flat energy time to calculate overall compressional energy.
TotalCompressionalEnergy = TotalEnergy(FlatEnergyIndex);
HSCompressionalEnergy = TotalHSEnergy(FlatEnergyIndex);
BangTimeCompressionEnergy = TotalEnergy(BangTimeIndex);
EnergyAtMaxPressure = TotalEnergy(CraxtonStagnationIndex);
EnergyAtMinHotspot = TotalEnergy(MaxOlsonCRIndex);


%Calculate estimated bang time alpha deposition as 1/2 of the overall alpha
%deposition, and use this to calculate the two f values.
%I looked at a few defns. Technically f should be bang time energy, but
%this doesnt work due to hotspot tracking. Q is stangnation, which also
%doesn't work. I am approximating Q as energy just after heating (although
%confusingly, I used the variable f at the time). Basically, ftot here is
%the value we want, and corresponds to Qtot in my electron heating paper.
EstimatedBangTimeAlphaDep = 0.5 * OverallTNEnergyDeposited;
ftot = EstimatedBangTimeAlphaDep/TotalCompressionalEnergy;
fhs = EstimatedBangTimeAlphaDep/HSCompressionalEnergy;
% ftot = EstimatedBangTimeAlphaDep/BangTimeCompressionEnergy;
% Qtot = EstimatedBangTimeAlphaDep/EnergyAtMinHotspot
end









if Plots==3;
      Intensity = (TotalLaserPower)/(4*pi()*(Radius(end,1)^2));
     figure
plot(Time./10^-9, Intensity./10^14, 'LineWidth', 2);
xlabel('Time (ns)')
ylabel('Intensity (10^{14} W/cm^2)')
%legend('Input Power', 'CBET Reduced Power', 'location', 'northwest')
set(gca,'yscale','log')
    
%Plot Pressure
figure
PressureGPa(PressureGPa<0) = 0;
PressureGPa(PressureGPa==0) = 0.1;
        surf(Time(1:end-1)*10^9,Radius(1:end-1,1:end-1), PressureGPa(:, 1:end-1));
        xlim([0 max(Time)]*10^9);
        ylim([0 max(Radius(:,1))]);
       set(gca,'ColorScale','log');
        colormap(jet);
        title('Pressure Log Plot');
        xlabel('Time (ns)');
        ylabel('Radius (cm)');
        zlabel('Pressure (GPa)');
        shading interp
        caxis([1 500])
        colorbar
        view(2)
        hold on
        plot3(Time*10^9,IceBoundary, max(PressureGPa), 'w');
        plot3(Time*10^9,CHBoundary, max(PressureGPa), 'w');
        hold off
        
%Plot zoning
figure
MassDifference = 100*diff(Mass(:,1))./Mass(2:end,1);
plot(MassDifference)
ylabel('Percentage Mass Difference')
yyaxis right;
plot(Mass(:,1), ':');
ylabel('Zone Mass')
xline(IceBoundaryIndex, ':');
xline(CHBoundaryIndex, ':');
title('Zoning Plot')
xlabel('Zone')

%Plot Laser Power
figure
plot(Time./10^-9, TotalLaserPower./10^12);
hold on
plot(Time./10^-9, SimulationPower./10^12, ':');
hold off;
title(['Total Energy = ', num2str(TotalLaserEnergy/1000), ' kJ, Required Energy = ', num2str(RequiredLaserEnergy/1000), ' kJ.'])
xlabel('Time (ns)')
ylabel('Laser Power (TW)')
legend('Required Power', 'Simulation Power')

MaxDensity = max(max(Density));

%Plot density, focussed around region of max density
figure;
DensityPlot = surf(Time,Radius(1:end-1,:), Density);
%xlim([(MaxDensityTime-.5e-9) (MaxDensityTime+.5e-9)]);
xlim([xlowerlim xupperlim]);
%ylim([0 0.01]);
ylim([0 yupperlim]);
zlim([0 MaxDensity]);
colorbar;
colormap(jet);
title('Density Plot')
xlabel('Time (ns)');
ylabel('Radius (cm)');
zlabel('Mass Density (g/cm^3)');
shading interp
view(2)
hold on
plot3(Time,IceBoundary, max(Density), 'w');
plot3(Time,CHBoundary, max(Density), 'w');
plot3(Time,HSRadius, max(Density),'linestyle', '--','color', 'w');
plot3(Time,ShellRadius, max(Density), 'linestyle', '--','color', 'w');
hold off

% %Plot OldRhoR, focussed around region of max density
% figure;
% OldRhoRPlot = surf(Time,Radius(1:end-1,:), OldRhoR);
% xlim([xlowerlim xupperlim]);
% ylim([0 yupperlim]);
% zlim([0 MaxOldRhoR]);
% colorbar;
% colormap(jet);
% title('OldRhoR Plot')
% xlabel('Time (ns)');
% ylabel('Radius (cm)');
% zlabel('OldRhoR (g/cm^2)');
% shading interp
% view(2)
% hold on
% plot3(Time,IceBoundary, max(OldRhoR), 'w');
% plot3(Time,CHBoundary, max(OldRhoR), 'w');
% hold off

%Plot Ion Temp, focussed around region of max density
figure
IonTempPlot = surf(Time,Radius(1:end-1,:), IonTemp);
xlim([xlowerlim xupperlim]);
ylim([0 yupperlim]);
%zlim([0 MaxIonTemp])
caxis([0 15000])
colorbar
colormap(jet)
title('Ion Temp Plot')
xlabel('Time (ns)')
ylabel('Radius (cm)')
zlabel('Ion Temperature (K)')
shading interp
view(2)
hold on
plot3(Time,IceBoundary, max(IonTemp), 'w');
plot3(Time,CHBoundary, max(IonTemp), 'w');
plot3(Time,HSRadius, max(IonTemp), 'linestyle', '--','color', 'w');
plot3(Time,ShellRadius, max(IonTemp), 'linestyle', '--','color', 'w');
hold off


%Simulate some data, and plot the Ion Temp vs OldRhoR against Cheng's analytic
%formula
figure
plot(HSAverageTemp(MinStagnationIndex:MaxStagnationIndex), HSAverageOldRhoR(MinStagnationIndex:MaxStagnationIndex));
hold on
plot(SimT, SimOldRhoR)
xlabel('Average Hot Spot Ion Temp (eV)')
ylabel('Average Hot Spot OldRhoR (g/cm^2)')
legend ('Simulated data', 'Cheng analytic formula')

% %Plot the average Rho R and average T for hotspot and shell around
% %stagnation
% figure
% tiledlayout(1,2);
% nexttile;
% title('Hotspot')
% yyaxis left;
% plot(Time./10^-9, HSAverageTemp);
% ylabel('Temperature, eV');
% xlim([Time(MinStagnationIndex)./10^-9 Time(MaxStagnationIndex)./10^-9])
% hold on;
% yyaxis right;
% plot(Time./10^-9, HSAverageOldRhoR);
% ylabel('OldRhoR, g/cm^2');
% xlim([Time(MinStagnationIndex)./10^-9 Time(MaxStagnationIndex)./10^-9])
% xline(BangTime./10^-9, 'k');
% xlabel('Time (ns)')
% nexttile;
% title('Shell')
% yyaxis left;
% plot(Time./10^-9, ShellAverageTemp);
% ylabel('Temperature, eV');
% xlim([Time(MinStagnationIndex)./10^-9 Time(MaxStagnationIndex)./10^-9])
% hold on;
% yyaxis right;
% plot(Time./10^-9, ShellAverageOldRhoR)
% ylabel('OldRhoR, g/cm^2');
% xlim([Time(MinStagnationIndex)./10^-9 Time(MaxStagnationIndex)./10^-9])
% xline(BangTime./10^-9, 'k');
% xlabel('Time (ns)')

if Type==1;
    %Plot TN reactions, focussed around region of max density
    figure
    IonTempPlot = surf(Time,Radius(1:end-1,:), TNOutput+0.001);
    xlim([xlowerlim xupperlim]);
    ylim([0 yupperlim]);
    set(gca,'ColorScale','log');
    colorbar
    colormap(jet)
    title('TN Reactions')
    xlabel('Time (ns)')
    ylabel('Radius (cm)')
    zlabel('TNOutput')
    shading interp
    view(2)
    hold on
    plot3(Time,IceBoundary, max(TNOutput), 'w');
    plot3(Time,CHBoundary, max(TNOutput), 'w');
    plot3(Time,HSRadius, max(TNOutput), 'linestyle', '--','color', 'w');
    plot3(Time,ShellRadius, max(TNOutput), 'linestyle', '--','color', 'w');
    hold off


    %Plot TN reactions, focussed around region of max density
    figure
    IonTempPlot = surf(Time,Radius(1:end-1,:), TNInput+0.001);
    xlim([xlowerlim xupperlim]);
    ylim([0 yupperlim]);
    set(gca,'ColorScale','log');
    colorbar
    colormap(jet)
    title('TN Deposition')
    xlabel('Time (ns)')
    ylabel('Radius (cm)')
    zlabel('TNInput')
    shading interp
    view(2)
    hold on
    plot3(Time,IceBoundary, max(TNInput), 'w');
    plot3(Time,CHBoundary, max(TNInput), 'w');
    plot3(Time,HSRadius, max(TNInput), 'linestyle', '--','color', 'w');
    plot3(Time,ShellRadius, max(TNInput), 'linestyle', '--','color', 'w');
    hold off
    
     
end

%Plot displaying shock trajectories (using inverse Pressure scale
     %length, as in Craxton's review paper
% InversePressureScaleLength = abs( diff(log(Pressure))./diff(Radius(1:end-1:end, :)));
% InversePressureScaleLength(isnan(InversePressureScaleLength)) = 0.001;
% InversePressureScaleLength(isinf(InversePressureScaleLength)) = 0.001; 
% 
%  figure
%         surf(Time*10^9,Radius(1:end-2,:), InversePressureScaleLength);
%         xlim([0 max(Time)*10^9]);
%         ylim([0 max(Radius(:,1))]);
%         set(gca,'ColorScale','log');
%         colormap(jet);
%         title('Density Log Plot');
%         xlabel('Time (ns)');
%         ylabel('Radius (cm)');
%         zlabel('Mass Density (g/cm^3)');
%         shading interp
%         colorbar
%                 view(2)
%         hold on
%         plot3(Time*10^9,IceBoundary, max(Density), 'k');
%         plot3(Time*10^9,CHBoundary, max(Density), 'k');
%         hold off
%         caxis([0.1 4])
% 
% % Plot inverse density scale length. Not sure on definition, but seems to
% % give a figure roughly like the one given by Tim Collins in emails. Need
% % to increase time res to make it look nice! Actually looks better if we
% don't multiply by density in first line.
InverseDensityScaleLength = abs( diff(log(Density))./diff(Radius(1:end-1:end, :)));
InverseDensityScaleLength(isnan(InverseDensityScaleLength)) = 0.001;
InverseDensityScaleLength((InverseDensityScaleLength==0)) = 0.001;
InverseDensityScaleLength(isinf(InverseDensityScaleLength)) = 0.001; 
Zones = repmat([1:size(Radius,1)-2].', 1, length(Time));

 figure
        surf(Time*10^9,Zones, InverseDensityScaleLength);
        xlim([0 max(Time)*10^9]);
        set(gca,'ColorScale','log');
        colormap(jet);
        title('Shock trajectories');
        xlabel('Time (ns)');
        ylabel('Cell number');
        zlabel('Mass Density (g/cm^3)');
        shading interp
        colorbar
                view(2)
        hold on
h = yline(CHBoundaryIndex, 'w', 'LineWidth', 2);
h = yline(IceBoundaryIndex, 'w', 'LineWidth', 2);

%Multiplot showing IFAR, Implosion Velocity, and Isentrope Parameter
%against time
figure
set(gcf, 'Position',  [100, 100, 2160, 520])
tiledlayout(1,3);
nexttile;
yyaxis left
plot(Time./10^-9, IFAR)
ylabel('IFAR');
yyaxis right
plot(Time./10^-9, abs(ImplosionVelocity))
ylabel('Velocity (km/s)');
xlabel('Time (ns)');
xline(BangTime./10^-9, 'k');
nexttile;
yyaxis left
plot(Time./10^-9, MinIsentrope)
ylabel('Min Isentrope Parameter');
ylim([0 9]);
yyaxis right
plot(Time./10^-9, AveragedIsentrope)
ylabel('Average Isentrope Parameter');
xlabel('Time (ns)');
ylim([0 9])
xline(BangTime./10^-9, 'k', {'Bang Time'});
xline(Time(MinStagnationIndex)./10^-9, '--k', {'0.001% Neutrons'})
nexttile;
yyaxis left
if isnan(max(Neutrons))~=1
plot(Time./10^-9, Neutrons);
xlim([0 max(Time)./10^-9]);
ylim([0 max(Neutrons)+1]);
xline(Time(i)./10^-9);
yline(Neutrons(i));
ylabel('Neutrons');
end
yyaxis right
plot(Time./10^-9, TotalLaserPower./10^12);
xlabel('Time (ns)');
ylabel('Laser Power (TW)');
xline(BangTime./10^-9, 'k');

end

if Plots==2 || Plots==3
    if sum(sum(isnan(Density)))==1
        %Log Density plot of overall full time and radius range
        figure
        surf(Time*10^9,Radius(1:end-1,:), Density);
        xlim([0 max(Time)*10^9]);
        ylim([0 max(Radius(:,1))]);
        set(gca,'ColorScale','log');
        colormap(jet);
        title('Density Log Plot');
        xlabel('Time (ns)');
        ylabel('Radius (cm)');
        zlabel('Mass Density (g/cm^3)');
        shading interp
        colorbar
        view(2)
        hold on
        plot3(Time*10^9,IceBoundary, max(Density), 'w');
        plot3(Time*10^9,CHBoundary, max(Density), 'w');
        hold off
    else
        figure
        surf(Time(1:end-1)*10^9,Radius(1:end-1,1:end-1), Density(:, 1:end-1));
        xlim([0 max(Time)]*10^9);
        ylim([0 max(Radius(:,1))]);
        set(gca,'ColorScale','log');
        colormap(jet);
        title('Density Log Plot');
        xlabel('Time (ns)');
        ylabel('Radius (cm)');
        zlabel('Mass Density (g/cm^3)');
        shading interp
        colorbar
        view(2)
        hold on
        plot3(Time*10^9,IceBoundary, max(Density), 'w');
        plot3(Time*10^9,CHBoundary, max(Density), 'w');
        hold off
    end

    %Density plot, but with zones displayed
    figure
    surf(Time,Radius(1:end-1,:), Density);
    xlim([0e-9 max(Time)]);
    ylim([0 0.08]);
    set(gca,'ColorScale','log');
    colormap(jet);
    title('Density Plot showing zone boundaries')
    xlabel('Time (ns)');
    ylabel('Radius (cm)');
    zlabel('Mass Density (g/cm^3)');
    colorbar
    view(2)
    hold on
    plot3(Time,IceBoundary, max(Density), 'w');
    plot3(Time,CHBoundary, max(Density), 'w');
    hold off
end

if Plots==1 || Plots==2 ||Plots==3
%Report max values/positions, and laser parameters
% disp(['Max power of ' ,num2str(max(TotalLaserPower)/10^12), ' TW, and total energy of ' num2str(TotalLaserEnergy/1000), ' kJ (Required: ' num2str(RequiredLaserEnergy/1000), ' kJ). Neutron energy of '  num2str(NeutronEnergy/1000) ' kJ. Gain of ' num2str(Gain)])
% disp(['Max density of ' ,num2str(MaxDensity), ' g/cm^3 at ', num2str(MaxDensityTime), ' s and a radius of ', num2str(MaxDensityRadius), ' cm (OldRhoR = ', num2str(MaxDensityOldRhoR), ', Ion Temp = ', num2str(MaxDensityTemp), ')'])
% disp(['Max OldRhoR of ' ,num2str(MaxOldRhoR), ' g/cm^2 at ', num2str(MaxOldRhoRTime), ' s and a radius of ', num2str(MaxOldRhoRRadius)])
% disp(['Max ion temperature of ' ,num2str(MaxIonTemp), ' eV at ', num2str(MaxIonTempTime), ' s and a radius of ', num2str(MaxIonTempRadius), ' cm(OldRhoR = ', num2str(MaxIonTempOldRhoR), ')'])

if max(OlsonCR)<20 || Type==0
    disp(['Olson CR of ',num2str(MaxOlsonCR) , ', Lindl CR of ',num2str(max(LindlCR(~isinf(LindlCR))))] )
else
    MeaningfulCRRange = find(diff(OlsonCR)>1,1, 'first');
    disp(['Olson CR of ',num2str(MaxOlsonCR) , ', Lindl CR of ',num2str(max(LindlCR(~isinf(LindlCR))))] )
    %[~, ShellMinIndex] = min(ShellRadius);
    %disp(['Convergence ratio of ' ,num2str(ConvergenceRatio), ', Olson CR is ROUGHLY ',num2str(OlsonCR(ShellMinIndex)) , ', Lindl CR is ROUGHLY ',num2str(LindlCR(ShellMinIndex))] )
end

disp(['IFAR of ' ,num2str(IFARCraxton)])
disp(['Implosion Velocity of ' ,num2str(abs(min(ImplosionVelocity)))])
fprintf('Parametric Limit of %e (Should be under 10^14) \n', ParametricLimit)
fprintf('%e neutrons produced \n', max(Neutrons))
disp(['Burn Up fraction of ' ,num2str(BurnUpFraction), '%'])
% disp(['Timesteps above ignition:' ,num2str(TimestepsAboveIgnition)])
disp(['Shell Kinetic Energy = ', num2str(ShellKineticEnergyValue/1000), ' kJ']);
disp(['Hydrodynamic Efficiency = ', num2str(HydroEfficiencyValue)]);

%Stagnation plot to show where CR was recorded
figure;
xlim([xlowerlim xupperlim]);
ylim([0 yupperlim]);
title('Stagnation Plot')
xlabel('Time (ns)');
ylabel('Radius (cm)');
hold on
plot(Time,IceBoundary, 'k');
plot(Time,CHBoundary,  'k');
plot(Time,HSRadius, 'linestyle', '--','color', 'k');
plot(Time,ShellRadius, 'linestyle', '--','color', 'k');
xline(BangTime, 'r');
xline(CraxtonStagnation, 'b');
plot(Time(MaxOlsonCRIndex),HSRadius(MaxOlsonCRIndex) ,'r*')
legend({'Ice Boundary', 'CH Boundary', 'HS Radius', 'Shell Radius', 'Bang Time', 'Stagnation'}, 'Location', 'NorthWest');
hold off

%Plot the average Rho R and average T for hotspot and shell around
%stagnation
figure
yyaxis left;
plot(Radius(1:end-1, BangTimeIndex)*10000, IonTemp(:, BangTimeIndex), 'LineWidth', 2);
ylabel('Temperature (eV)');
hold on;
yyaxis right;
plot(Radius(1:end-1, BangTimeIndex)*10000, Density(:, BangTimeIndex), 'LineWidth', 2);
ylabel('Density (g/cm^3)');
xlabel('Radius (µm)')
xlim([0 ShellRadius(BangTimeIndex)*1.5*10000])
xline(HSRadius(BangTimeIndex)*10000,'linestyle', '--','color', 'k', 'LineWidth', 2);
xline(ShellRadius(BangTimeIndex)*10000, 'linestyle', '--','color', 'k', 'LineWidth', 2);
title('At bang time')
fig=gcf;
fig.Units               = 'points';
fig.Position(3)         = 300;
fig.Position(4)         = 200;
set(fig.Children, ...
    'FontName',     'Times', ...
    'FontSize',     9);
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))


%Plot the average Rho R and average T for hotspot and shell around
%stagnation
figure
yyaxis left;
plot(Radius(1:end-1, CraxtonStagnationIndex)*10000, IonTemp(:, CraxtonStagnationIndex), 'LineWidth', 2);
ylabel('Temperature (eV)');
hold on;
yyaxis right;
plot(Radius(1:end-1, CraxtonStagnationIndex)*10000, Density(:, CraxtonStagnationIndex), 'LineWidth', 2);
ylabel('Density (g/cm^3)');
xlabel('Radius (µm)')
xlim([0 ShellRadius(CraxtonStagnationIndex)*1.5*10000])
xline(HSRadius(CraxtonStagnationIndex)*10000,'linestyle', '--','color', 'k', 'LineWidth', 2);
xline(ShellRadius(CraxtonStagnationIndex)*10000, 'linestyle', '--','color', 'k', 'LineWidth', 2);
title('At stagnation')
fig=gcf;
fig.Units               = 'points';
fig.Position(3)         = 300;
fig.Position(4)         = 200;
set(fig.Children, ...
    'FontName',     'Times', ...
    'FontSize',     9);
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))

if isempty(HydroEfficiencyValue)~=1 && isnan(HydroEfficiencyValue)~=1
figure
plot(Time, HydroEfficiency);
hold on
xline(BangTime, '--r');
xline(CraxtonStagnation, '--b');
xline(Time(find(ImplosionVelocity==min(ImplosionVelocity), 1, 'first')), '--k');
yline(HydroEfficiencyValue);
hold off
legend({'Efficiency', 'Bang Time', 'Stagnation', 'Max Velocity'}, 'Location', 'NorthWest');
title('Hydrodynamic Efficiency')
ylim([0 1.5*HydroEfficiencyValue])
end

Report_Git

figure
tiledlayout(1,2);
nexttile;
title('Hotspot')
yyaxis left;
plot(Time./10^-9, HSAverageTemp);
ylabel('Temperature, eV');
xlim([Time(MinStagnationIndex)./10^-9 Time(MaxStagnationIndex)./10^-9])
hold on;
yyaxis right;
plot(Time./10^-9, HSRhoR);
ylabel('RhoR, g/cm^2');
xlim([Time(MinStagnationIndex)./10^-9 Time(MaxStagnationIndex)./10^-9])
xline(BangTime./10^-9, 'k');
xlabel('Time (ns)')
nexttile;
title('Shell')
yyaxis left;
plot(Time./10^-9, ShellAverageTemp);
ylabel('Temperature, eV');
xlim([Time(MinStagnationIndex)./10^-9 Time(MaxStagnationIndex)./10^-9])
hold on;
% yyaxis right;
% plot(Time./10^-9, CapsuleRhoR)
% hold on
% plot(Time./10^-9, CapsuleRhoRDTOnly, '--')
% hold off
% ylabel('Overall Capsule RhoR, g/cm^2');
% xlim([Time(MinStagnationIndex)./10^-9 Time(MaxStagnationIndex)./10^-9])
% xline(BangTime./10^-9, 'k');
% xlabel('Time (ns)')
yyaxis right;
plot(Time./10^-9, ShellRhoR)
ylabel('\rho R, g/cm^2');
xlim([Time(MinStagnationIndex)./10^-9 Time(MaxStagnationIndex)./10^-9])
xline(BangTime./10^-9, 'k');
xlabel('Time (ns)')

try
  figure
plot(Time.*10^9, TotalThermalEnergy./1000);
hold on
plot(Time.*10^9, TotalKineticEnergy./1000);
plot(Time.*10^9, TotalEnergy./1000, 'k');
xline(BangTime*10^9, '--k');
xline(Time(MaxVelocityIndex)*10^9, '--k');
xline(Time(FlatEnergyIndex)*10^9, '--r');
end
end

