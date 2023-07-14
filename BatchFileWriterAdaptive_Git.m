%Script to write Hyades input decks. Can write for three or four layers
%(specified by Regions) and three or 4 pulses.


%Output timing
%ChangeTimeAt = [0 16 22]*10^-9;
%TimeSpacing = [1E-10 1E-11];
Time=[0];
for i=1:(length(ChangeTimeAt)-1)
    Time = [Time, (ChangeTimeAt(i)+TimeSpacing(i)):TimeSpacing(i):ChangeTimeAt(i+1)];
end


%Densities of the three materials
DensityIce = 2.530000e-01; %0.212 in Helios, 0.253 in Robbie's
DensityCH = 1.040000e+00; %1.11 in Helios, 1.04 in Robbie's
% DensityCH = 2.530000e-01; %0.212 in Helios, 0.253 in Robbie's

IonTempAll = 1.551000e-06;

%Calculates relevant boundaries from layers
if Regions==3
BoundaryVapour = LayersVapour+1;
BoundaryIce = BoundaryVapour + LayersIce;
BoundaryCH = BoundaryIce + LayersCH;
else
BoundaryVapour = LayersVapour+1;
BoundaryIce = BoundaryVapour + LayersIce;
BoundaryCH1 = BoundaryIce + LayersCH1;
BoundaryCHEndZone = BoundaryCH1 +LayersCHEndZone;
end

%Converts Pulse Power in TW into erg, and reduces to account for cross beam
%energy transfer
Pulse1erg = Pulse1Power * 0.8*10^19;
Pulse2erg = Pulse2Power * 0.8*10^19;
Pulse3erg = Pulse3Power * 0.8*10^19;
if exist('Pulse4Time','var')==1
Pulse4erg = Pulse4Power * 0.8*10^19;
end
if SecondLaser==1
SecondLasererg = SecondLaserPower * 0.8*10^19;
end

%Print the file
filename = sprintf(['File', num2str(ArrayIndex)]);
filepath = sprintf(['MultiFile/', filename, '.inf']);
fileID = fopen(filepath,'w');

fprintf(fileID,'%s\r\nc\r\n', filename);
fprintf(fileID,'geometry 3 \r\nc\r\n');

if Regions==3
fprintf(fileID,'mesh 1 %.0f 0 %.6e %.6e \r\n', BoundaryVapour, RadiusVapour, RatioVapour);
fprintf(fileID,'mesh %.0f %.0f %.6e %.6e %.6e \r\n', BoundaryVapour, BoundaryIce, RadiusVapour, RadiusIce, RatioIce);
fprintf(fileID,'mesh %.0f %.0f %.6e %.6e %.6e \r\nc\r\n', BoundaryIce, BoundaryCH, RadiusIce, RadiusCH, RatioCH);
else
fprintf(fileID,'mesh 1 %.0f 0 %.6e %.6e \r\n', BoundaryVapour, RadiusVapour, RatioVapour);
fprintf(fileID,'mesh %.0f %.0f %.6e %.6e %.6e \r\n', BoundaryVapour, BoundaryIce, RadiusVapour, RadiusIce, RatioIce);
fprintf(fileID,'mesh %.0f %.0f %.6e %.6e %.6e \r\n', BoundaryIce, BoundaryCH1, RadiusIce, RadiusCHEndZone, RatioCH1);
fprintf(fileID,'mesh %.0f %.0f %.6e %.6e %.6e \r\nc\r\n', BoundaryCH1, BoundaryCHEndZone, RadiusCHEndZone, RadiusCH, RatioCHEndZone);
end

if Regions==3
fprintf(fileID, 'region 1 %.0f 1 %.3e %.3e \r\n', LayersVapour, DensityVapour, IonTempAll);
fprintf(fileID, 'region %.0f %.0f 2 %.3e %.3e \r\n', BoundaryVapour, (BoundaryIce-1), DensityIce, IonTempAll);
fprintf(fileID, 'region %.0f %.0f 3 %.3e %.3e \r\nc\r\n', BoundaryIce, (BoundaryCH-1), DensityCH, IonTempAll);
else
fprintf(fileID, 'region 1 %.0f 1 %.3e %.3e \r\n', LayersVapour, DensityVapour, IonTempAll);
fprintf(fileID, 'region %.0f %.0f 2 %.3e %.3e \r\n', BoundaryVapour, (BoundaryIce-1), DensityIce, IonTempAll);
fprintf(fileID, 'region %.0f %.0f 3 %.3e %.3e \r\n', BoundaryIce, (BoundaryCH1-1), DensityCH, IonTempAll);
fprintf(fileID, 'region %.0f %.0f 4 %.3e %.3e \r\nc\r\n', BoundaryCH1, (BoundaryCHEndZone-1), DensityCH, IonTempAll);
end

fprintf(fileID, 'material 1 dt \r\n');
if FoamOnly==1
    fprintf(fileID, 'material 2 1. 2.014  0.5 \r\n');
fprintf(fileID, 'material 2 6. 12.012 0.5 \r\n');
else 
    fprintf(fileID, 'material 2 dt \r\n');
end
if Regions==3
fprintf(fileID, 'material 3 1. 2.014  0.5 \r\n');
fprintf(fileID, 'material 3 6. 12.012 0.5 \r\nc\r\n');
else
    fprintf(fileID, 'material 3 1. 2.014  0.5 \r\n');
fprintf(fileID, 'material 3 6. 12.012 0.5 \r\n');
fprintf(fileID, 'material 4 1. 2.014  0.5 \r\n');
fprintf(fileID, 'material 4 6. 12.012 0.5 \r\nc\r\n');
end

fprintf(fileID, 'qeos 1 1.86e9 3.0e-04 \r\n');
if FoamOnly==1
    fprintf(fileID, 'eos /work3/clf/rad_hydro/hyades/EOS-Opacity/SESAME//eos_32.dat 2 \r\n');
else
    fprintf(fileID, 'qeos 2 1.86e9 0.253 \r\n');
end
if Regions==3
fprintf(fileID, 'eos /work3/clf/rad_hydro/hyades/EOS-Opacity/SESAME//eos_32.dat 3 \r\nc\r\n');
else
    fprintf(fileID, 'eos /work3/clf/rad_hydro/hyades/EOS-Opacity/SESAME//eos_32.dat 3 \r\n');
    fprintf(fileID, 'eos /work3/clf/rad_hydro/hyades/EOS-Opacity/SESAME//eos_32.dat 4 \r\nc\r\n');
end

fprintf(fileID, 'eosxtrp  1  1  2  1  2 \r\n');
fprintf(fileID, 'eosxtrp  2  1  2  1  2 \r\n');
if Regions==3
fprintf(fileID, 'eosxtrp  3  1  2  1  2 \r\nc\r\n');
else
    fprintf(fileID, 'eosxtrp  3  1  2  1  2 \r\n');
fprintf(fileID, 'eosxtrp  4  1  2  1  2 \r\nc\r\n');
end

fprintf(fileID, 'ioniz 1 4 \r\n');
fprintf(fileID, 'ioniz 2 4 \r\n');
if Regions==3
fprintf(fileID, 'ioniz 3 4 \r\nc\r\n');
else
    fprintf(fileID, 'ioniz 3 4 \r\n');
fprintf(fileID, 'ioniz 4 4 \r\nc\r\n');
end

if Regions==3
fprintf(fileID, 'source laser %.3f -%.0f \r\n', Wavelength, BoundaryCH);

else 
    fprintf(fileID, 'source laser %.3f -%.0f \r\n', Wavelength, BoundaryCHEndZone);

end
fprintf(fileID, 'tv 0. 0. \r\n');
fprintf(fileID, 'tv %.3e %.3e \r\n', RiseTime, Pulse1erg);
fprintf(fileID, 'tv %.3e %.3e \r\n', Pulse2Time, Pulse1erg);
fprintf(fileID, 'tv %.3e %.3e \r\n', (Pulse2Time+RiseTime), Pulse2erg);
fprintf(fileID, 'tv %.3e %.3e \r\n', Pulse3Time, Pulse2erg);
fprintf(fileID, 'tv %.3e %.3e \r\n', (Pulse3Time+RiseTime), Pulse3erg);
if exist('Pulse4Time','var')==0
fprintf(fileID, 'tv %.3e %.3e \r\n', (LaserEndTime), Pulse3erg);
else
    fprintf(fileID, 'tv %.3e %.3e \r\n', Pulse4Time, Pulse3erg);
    fprintf(fileID, 'tv %.3e %.3e \r\n', (Pulse4Time+RiseTime), Pulse4erg);
    fprintf(fileID, 'tv %.3e %.3e \r\n', (LaserEndTime), Pulse4erg);
end
fprintf(fileID, 'tv %.3e 0. \r\nc\r\n', (LaserEndTime+RiseTime));



if SecondLaser==1
if Regions==3
fprintf(fileID, 'source laser %.3e -%.0f \r\n', SecondLaserWavelength, BoundaryCH);
else 
    fprintf(fileID, 'source laser %.3e -%.0f \r\n', SecondLaserWavelength, BoundaryCHEndZone);
end
fprintf(fileID, 'tv 0. 0. \r\n');
fprintf(fileID, 'tv %.3e 0. \r\n', SecondLaserOnTime);
fprintf(fileID, 'tv %.3e %.3e \r\n', (SecondLaserOnTime+RiseTime), SecondLasererg);
fprintf(fileID, 'tv %.3e %.3e \r\n', SecondLaserOffTime, SecondLasererg);
fprintf(fileID, 'tv %.3e 0. \r\nc\r\n', (SecondLaserOffTime+RiseTime));
end




if ElectronHeatingOn == 1
    fprintf(fileID, 'source ee 1 %.0f \r\n', HeatingZones-1);
    fprintf(fileID, 'tv 0. 0. \r\n');
    for i=1:1:length(HeatingTimings)
        fprintf(fileID, 'tv %.5e %.5e \r\n', HeatingTimings(i), Heating(i));
    end
end
fprintf(fileID, 'c\r\n');



if TNBurn == 1
    if FoamOnly==1
    fprintf(fileID, 'tnburn 1 \r\n');
    else
    fprintf(fileID, 'tnburn 1 2 \r\n');
    end
    fprintf(fileID, 'parm tibmn 0.1  \r\nc\r\n');
else
    fprintf(fileID, 'c No TNBurn in this run \r\nc\r\n');
end

fprintf(fileID, 'group  0 20 0.03 1.0 \r\n');
fprintf(fileID, 'group 20 50 1.00 5.0 \r\n');
fprintf(fileID, 'group 50 70 5.00 300.0 \r\nc\r\n');

fprintf(fileID, 'pparray rho te ti tr pres R Rcm zbar u deplas xlsint vol bpeprd bpeprdr bpedep dene eion eelc TNDENI \r\nc\r\n');
fprintf(fileID, 'parm xlibam 1.0 \r\n');
fprintf(fileID, 'parm flxlem 0.050 \r\n');
fprintf(fileID, 'parm flxlim 0.4 \r\n');
fprintf(fileID, 'parm alvism 0.3 \r\n');
fprintf(fileID, 'parm aqvism 2.0 \r\n');
fprintf(fileID, 'parm qstimx 4.3e-6 \r\n');
fprintf(fileID, 'parm lrdtrn 1 \r\n');
fprintf(fileID, 'parm temin 1.551e-06 \r\n');
fprintf(fileID, 'parm timin 1.551e-06 \r\n');
fprintf(fileID, 'parm irdtrn 2 \r\n');
fprintf(fileID, 'parm nstop 1e8 \r\n');
fprintf(fileID, 'parm dt 1e-15 \r\n');
fprintf(fileID, 'parm dtmin 1e-25  \r\n');
fprintf(fileID, 'parm JHTRMX 200 \r\n');
for i=2:1:(length(ChangeTimeAt)-1)
    fprintf(fileID, 'change %.4e postdt %.4e \n', ChangeTimeAt(i), TimeSpacing(i));
end
fprintf(fileID,'parm postdt %.4e \n', TimeSpacing(1));
fprintf(fileID,'parm tstop %.4e \n', ChangeTimeAt(end));




fclose(fileID);





%Function to recreate the Hyades mesh function - given the mesh boundaries
%j and radii to stretch between r and the ratio, it will create the mesh
%coordinates.
function Radii = RatioIncrement(j1, j2, r1, r2, Ratio)
Index = (1+j1-j1):(j2-j1);
Increment = Ratio.^(Index);
Radius = cumsum(Increment);
Scaling = (r2-r1)/Radius(end);
Radii = [r1+(Radius*Scaling)];
end

%Given the radii and density, will calculate the mass difference between
%zones and return the mean of the square of this value.
function [MaxDiff, ZoneMass, MassDiff] = MassDifference(Radii, Density)
Vol = (4/3) * pi() * (Radii.^3);
ZoneMass = diff(Vol).*Density;
MassDiff = 100*(diff(ZoneMass)./ZoneMass(2:end));
%Good results for the mean of the quadrature. This function takes the
%mean of quadrature of only those values with a mass difference above
%1.5 (bear in mind the vapour layer is not optimised, so early high
%values are constant).
MaxDiff = mean((abs(MassDiff).*(abs(MassDiff)>2)).^2);
end

%Optimisation function. Given the Ratios, Radii and Layers, will construct
%the mesh and calculate the mass difference using above functions.
function [MaxDiff, Radii, Density] = BoundarySolver(RatioVapour, RatioSplit, RatioIce, RatioCH, LayersVapour, LayersSplit, LayersIce, LayersCH, RadiusVapour,RadiusSplit, RadiusIce, RadiusCH, DensityVapour, DensityIce, DensityCH)

BoundaryVapour = LayersVapour+1;
BoundarySplit = BoundaryVapour + LayersSplit;
BoundaryIce = BoundarySplit + LayersIce;
BoundaryCH = BoundaryIce + LayersCH;

VapourRadii = RatioIncrement(1, BoundaryVapour, 0, RadiusVapour, RatioVapour);
SplitRadii = RatioIncrement(BoundaryVapour, BoundarySplit, RadiusVapour, RadiusSplit, RatioSplit);
IceRadii = RatioIncrement(BoundarySplit, BoundaryIce, RadiusSplit, RadiusIce, RatioIce);
CHRadii = RatioIncrement(BoundaryIce, BoundaryCH, RadiusIce, RadiusCH, RatioCH);

Density = [DensityVapour*ones(1, LayersVapour), DensityIce*ones(1, (LayersIce+LayersSplit)), DensityCH*ones(1, LayersCH)];

Radii = [0, VapourRadii, SplitRadii, IceRadii, CHRadii];
MaxDiff = MassDifference(Radii, Density);

end


