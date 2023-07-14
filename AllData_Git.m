%Choose mode. 1 plots for all data, while 2 allows a criteria to be
%specified. 0 gives all the stored information on a given data point
Mode=2;

%Choose which data sets to include
Include3Pulse=0;
Include4Pulse=1;
IncludeFoamOnly=0;
IncludeSecondLaser=0;

%Compare3and4=1 makes these different colours, while off treats them the
%same and uses colour for different variables.
Compare3and4=1;

Paths=[];
Files=[];
%Set Paths for Data where MatlabWorkspace is saved
if Include3Pulse==1;
Import = dir('\\alfs1.physics.ox.ac.uk\al\paddock\Hyades Data\*.mat');
Paths3Pulse =  string({Import.name}.');
Paths = [Paths; Paths3Pulse];
Files = [Files; strcat('\\alfs1.physics.ox.ac.uk\al\paddock\Hyades Data\', Paths3Pulse)];
end
if Include4Pulse==1;
Import = dir('\\alfs1.physics.ox.ac.uk\al\paddock\Hyades Data\4 Pulses\*.mat');
Paths4Pulse =  string({Import.name}.');
Paths = [Paths; Paths4Pulse];
Files = [Files; strcat('\\alfs1.physics.ox.ac.uk\al\paddock\Hyades Data\4 Pulses\', Paths4Pulse)];
end
if IncludeFoamOnly==1;
Paths=[];
Files=[];
Import = dir('\\alfs1.physics.ox.ac.uk\al\paddock\Hyades Data\Foam Only\*.mat');
PathsFoamOnly =  string({Import.name}.');
Paths = [Paths; PathsFoamOnly];
Files = [Files; strcat('\\alfs1.physics.ox.ac.uk\al\paddock\Hyades Data\Foam Only\', PathsFoamOnly)];
end
if IncludeSecondLaser==1;
Paths=[];
Files=[];
Import = dir('\\alfs1.physics.ox.ac.uk\al\paddock\Hyades Data\Second Laser\*.mat');
PathsSecondLaser =  string({Import.name}.');
Paths = [Paths; PathsSecondLaser];
Files = [Files; strcat('\\alfs1.physics.ox.ac.uk\al\paddock\Hyades Data\Second Laser\', PathsSecondLaser)];
end


%Load in variables, and in to one large array over all data
load(Files(1),'RequiredLaserEnergyAll','GainAll', 'OlsonCRAll', 'ValidFiles', ...
     'HydroEfficiencyAll', 'TotalEfficiencyAll', 'VapourRadiusAll', 'IceRadiusAll', 'CHRadiusAll', ...
     'VapourPressureAll', 'PulseAll', 'PowerAll', 'ShellKineticEnergyAll', 'IFARAll', 'ImplosionVelocityAll', 'NeutronsAll')
 if size(PulseAll,1)==3
     PulseAll = [PulseAll; zeros(1, size(PulseAll,2))];
PowerAll = [PowerAll; zeros(1, size(PowerAll,2))];
 end
 
OverallEnergy = RequiredLaserEnergyAll;
OverallGain = GainAll;
OverallCR = OlsonCRAll;
OverallImplosionVelocity = ImplosionVelocityAll;
OverallIFAR = IFARAll;
OverallValidFiles = ValidFiles;
OverallPulseTimes = PulseAll;
OverallPulsePower = PowerAll;
OverallShellKinetic = ShellKineticEnergyAll;
OverallHydroEfficiency = HydroEfficiencyAll;
OverallTotalEfficiency = TotalEfficiencyAll;
OverallVapourRadius = VapourRadiusAll;
OverallIceRadius = IceRadiusAll;
OverallCHRadius = CHRadiusAll;
OverallVapourPressure = VapourPressureAll;
OverallNeutrons = NeutronsAll;

OverallPath = repmat(Paths(1),1, length(GainAll));
DataSeries = repmat(1,1, length(GainAll));

for i=2:length(Files)
    i;
    load(Files(i),'RequiredLaserEnergyAll','GainAll', 'OlsonCRAll', 'ValidFiles', ...
     'HydroEfficiencyAll', 'TotalEfficiencyAll', 'VapourRadiusAll', 'IceRadiusAll', 'CHRadiusAll', ...
     'VapourPressureAll', 'PulseAll', 'PowerAll', 'ShellKineticEnergyAll', 'IFARAll', 'ImplosionVelocityAll', 'NeutronsAll')
  if size(PulseAll,1)==3
     PulseAll = [PulseAll; zeros(1, size(PulseAll,2))];
PowerAll = [PowerAll; zeros(1, size(PowerAll,2))];
 end
    OverallEnergy = [OverallEnergy, RequiredLaserEnergyAll];
    OverallGain = [OverallGain, GainAll];
    OverallNeutrons = [OverallNeutrons, NeutronsAll];
    OverallCR = [OverallCR, OlsonCRAll];
    OverallImplosionVelocity = [OverallImplosionVelocity, ImplosionVelocityAll];
    OverallIFAR = [OverallIFAR, IFARAll];
    OverallValidFiles = [OverallValidFiles, ValidFiles];
    OverallPulseTimes = [OverallPulseTimes, PulseAll];
    OverallPulsePower = [OverallPulsePower, PowerAll];
    OverallShellKinetic = [OverallShellKinetic,ShellKineticEnergyAll];
    OverallHydroEfficiency = [OverallHydroEfficiency,HydroEfficiencyAll];
    OverallTotalEfficiency = [OverallTotalEfficiency,TotalEfficiencyAll];
    OverallVapourRadius = [OverallVapourRadius,VapourRadiusAll];
    OverallIceRadius = [OverallIceRadius,IceRadiusAll];
    OverallCHRadius = [OverallCHRadius,CHRadiusAll];
    OverallVapourPressure = [OverallVapourPressure, VapourPressureAll];
    OverallPath = [OverallPath, repmat(Paths(i),1, length(GainAll))];
    DataSeries = [DataSeries, repmat(i,1, length(GainAll))];
end

%For All data plot, plot Gain against Energy, with CR giving the colour
if Mode==1
fig1 = figure('DeleteFcn','datacursormode');
    pointsize = 20;
    title('Neutrons (consider end time!)')
    scatter(OverallEnergy./(10^6),OverallGain,pointsize, OverallCR);
    xlabel('Input Energy (MJ)');
    ylabel('Gain')
    colorbar()
    dcm_obj = datacursormode(fig1);
    set(dcm_obj,'UpdateFcn',{@myupdatefcn,OverallCR, OverallPath, DataSeries, OverallValidFiles});
    colormap(copper)
    caxis([13 17])
end

%Set a criteria. New variables are created, containing the data for just
%those points that fit the criteria. Plot.
if Mode==2
    Criteria1 = OverallCR<16;
    Criteria2 = OverallIFAR<30;
    Criteria3 = OverallImplosionVelocity<400;
    Criteria = logical(Criteria1 .* Criteria2 .* Criteria3);
    %Criteria = logical(Criteria1 .* Criteria2);
    %Criteria = ((DataSeries == 12) + (DataSeries == 13) +(DataSeries == 14))>0.5;
    EnergyPlot = OverallEnergy(Criteria);
    GainPlot = OverallGain(Criteria);
    CRPlot = OverallCR(Criteria);
    ValidFilesPlot = OverallValidFiles(Criteria);
    PathPlot = OverallPath(Criteria);
    DataSeriesPlot = DataSeries(Criteria);
    fig1 = figure('DeleteFcn','datacursormode');
    pointsize = 20;
    if Compare3and4==1
        ColourCompare = OverallPulsePower(4,Criteria)>0;
        scatter(EnergyPlot./(10^6), GainPlot, pointsize, ColourCompare);
        colormap([1 0 0; 0 0 1]);
        caxis([0 1]);    
    else
    scatter(EnergyPlot./(10^6), GainPlot, pointsize, CRPlot);
    colormap(copper)
    caxis([14 16])
    end
    xlabel('Input Energy  (MJ)');
    ylabel('Gain')
    colorbar()
    dcm_obj = datacursormode(fig1);
    set(dcm_obj,'UpdateFcn',{@myupdatefcn,CRPlot, PathPlot, DataSeriesPlot, ValidFilesPlot});
  end
    
 if Mode == 0
     PathIdentifier = 24;
     FileNumber = 102;
     
     Criteria = (DataSeries == PathIdentifier).*(OverallValidFiles == FileNumber);
     FileIndex = find(Criteria == 1, 1);
     disp(strcat('File =  ', OverallPath(FileIndex), '\File', num2str(OverallValidFiles(FileIndex))))
     disp(strcat('Gain =  ', num2str(OverallGain(FileIndex))))
     disp(strcat('Energy =  ', num2str(OverallEnergy(FileIndex)/1000000), ' MJ'))
     disp(strcat('CR =  ', num2str(OverallCR(FileIndex))))
     disp(strcat('IFAR =  ', num2str(OverallIFAR(FileIndex))))
     disp(strcat('Implosion Velocity =  ', num2str(OverallImplosionVelocity(FileIndex)), ' km/s'))
     disp(strcat('HydroEfficiency =  ', num2str(OverallHydroEfficiency(FileIndex)), '%, TotalEfficiency = ', num2str(OverallTotalEfficiency(FileIndex)), '%, Shell KE = ', num2str(OverallShellKinetic(FileIndex)/1000), ' kJ'))
     disp(['Vapour Radius = ', num2str(OverallVapourRadius(FileIndex))]);
     disp(['Ice Radius = ', num2str(OverallIceRadius(FileIndex))]);
    disp(['CH Radius = ', num2str(OverallCHRadius(FileIndex))]);
    disp(['Vapour Pressure = ', num2str(OverallVapourPressure(FileIndex))]);
     disp(strcat( 'Pulse Times = ', num2str((abs(OverallPulseTimes(:, FileIndex)./10^-9)).', '%f ns,  '), '    Pulse Powers = ', num2str((abs(OverallPulsePower(:, FileIndex))).', '%f TW,  ')))

 end
 
    function txt = myupdatefcn(~,event_obj,OverallCR, OverallPath, DataSeries, OverallValidFiles)
% Customizes text of data tips
pos = get(event_obj,'Position');
I = get(event_obj, 'DataIndex');
txt = {['X: ',num2str(pos(1))],...
       ['Y: ',num2str(pos(2))],...
       ['CR: ',num2str(OverallCR(I))],...
       [strcat(OverallPath(I),'\File',num2str(OverallValidFiles(I)))],...
       ['Path Identifier: ', num2str(DataSeries(I))]};
end