  %Meshing script. Attempts to mesh using 4 layers, with a split in the CH
%region. If this fails or the split occurs too close to the Ice/CH
%interface, then it reverts to a layer mesh. The split is chosen so that
%the percentage mass difference of the CH end zone equals that of the ice
%layer. Changes now made for higher pressures

%Figures==1 sets manual mode - specify the capsule dimensions below, and
%figures will be plotted for the meshing.
Figures =0;

if Figures==1
    close all
    clear all
    
    Figures = 1;
    
   %Select Values for Radii to vary over
    RadiusVapour = 0.0895;
    RadiusIce = 0.0976;
    RadiusCH= 0.09975;
    DensityVapour=[0.00135];
    
end

%Determine size relative to a default value (this will be used in a scaling
%later to maintain accuracy regardless of size).
SizeMod = RadiusCH/0.24225;

%Insert Densities of the materials (if not specified in another script).
DensityIce = 2.530000e-01;
DensityCH = 1.040000e+00;

%Declare any constant thicknesses at the boundary. Vacuum thickness should
%be 0.01um (i.e. thickness of ch at vacuum)
ThicknessAtVacuum=0.000001; % 0.000001 Corresponds to 0.01um
%  ThicknessAtVacuum=0.00000001;

%Declare ratio and layers for DT vapour
LayersVapour = 150; %150
RatioVapour=1;

%Calculate Vapour radii and masses, and find the mass of the last vapour
%zone. This mass will be used as the target mass for our first ice zone.
BoundaryVapour = LayersVapour+1;
VapourRadii = RatioIncrement(1, BoundaryVapour, 0, RadiusVapour, RatioVapour);
VapourRadii = [0, VapourRadii];
VapourDensity = DensityVapour*ones(1, LayersVapour);
[VapourZoneMass, VapourMassDiff] = MassDifference(VapourRadii, VapourDensity);
VapourEndMass = VapourZoneMass(end);
ThicknessAtVapour = VapourEndMass/(4*pi()*DensityIce*(RadiusVapour^2));

%Calculate the thickness for the ice/CH interface. Determined a value for
%this which gives the desired percentage mass difference, and multiply by
%various parameters (determined through reasoning and trial and error) to
%attempt to remove their influence and make scale invariant.
%A smaller values gives a lower percentage mass difference, but required a
%larger number of layers. We find still that as the ice region grows the
%percentage mass may need to be reduced to avoid errors.
%Calculate this value using both the CH region and the ice region. As we
%want to maintain a maximum percentage error, we choose the smaller of
%these two values as the thickness to use.
%I believe that the two 4 layer systems should only ever use the Ice
%version of this value. The CH values may also need sizemod and vapour
%density multipliers, but I have not experimented with this yet.
%Whatever thickness is chosen for the ice, this must be multiplied by the
%density modifiers to give a thickness for CH that gives the same zone
%mass.
ThicknessAtCHForIceUsingIce =4.5e+03 *(1/SizeMod)* (RadiusIce-RadiusVapour) * ThicknessAtVapour *1/(DensityVapour/0.00065);
ThicknessAtCHForIceUsingCH = 90e+03 * (RadiusCH-RadiusIce) * 0.000001 / (DensityIce/DensityCH) *(1/SizeMod) ; %Was*10? %*1/(DensityVapour/0.00065)
ThicknessAtCHForIce = min([ThicknessAtCHForIceUsingIce, ThicknessAtCHForIceUsingCH]);
ThicknessAtCHForCH = ThicknessAtCHForIce*(DensityIce/DensityCH);

%We now try in sequence 3 approaches.
%1) The first is a four layer approach with a split in the CH layer. The
%ice-CH thickness is chosen above. The thickness at the split is the same
%as this thickness. The code finds the linear gradient in thickness vs.
%radius for the ice region, and using this estimates the radius range
%required for the CH to get from the vacuum thickness at the end to this
%thickness, maintaining the same percentage mass difference as the ice
%layer. The intermediate CH zone then maintains the same (large) thickness,
%reducing the number of zones. If this fails, we try the next approach.
try
    ThicknessAtCHForIceUsingIce =4.5e+03*(RadiusIce-RadiusVapour)*(1/SizeMod) * ThicknessAtVapour *1/(DensityVapour/0.00065); %4.5e+03
    ThicknessAtCHForIce = min([ThicknessAtCHForIceUsingIce, ThicknessAtCHForIceUsingCH]);
    ThicknessAtCHForCH = ThicknessAtCHForIce*(DensityIce/DensityCH);
    
    %Solve for Ice Layer
    [RatioIce, LayersIce] = ThicknessSolver(RadiusVapour, RadiusIce, ThicknessAtVapour, ThicknessAtCHForIce);
    
    %Calculate thickness for CH end zone, and gradient from the ice region.
    %Use this to calculate radius that the end zone starts at
    ThicknessAtCHEndZone = ThicknessAtCHForCH*(RadiusIce^2)/(RadiusCH^2);
    Gradient = (ThicknessAtCHForIce - ThicknessAtVapour)/(RadiusIce - RadiusVapour);
    RadiusCHEndZone = RadiusCH - ((ThicknessAtCHEndZone - ThicknessAtVacuum)/Gradient);
    RadiusCHEndZone = round(RadiusCHEndZone,4);
    
    
    %Solve for the first CH zone
    [RatioCH1, LayersCH1] = ThicknessSolver(RadiusIce, RadiusCHEndZone, ThicknessAtCHForCH, ThicknessAtCHEndZone);
    %Solve for the end zone
    [RatioCHEndZone, LayersCHEndZone] = ThicknessSolver(RadiusCHEndZone, RadiusCH, ThicknessAtCHEndZone, ThicknessAtVacuum);
    
    %Combine all and plot
    [Radii, Density] = BoundarySolver4Layer(RatioVapour, round(RatioIce, 4), round(RatioCH1, 4),  round(RatioCHEndZone, 4),LayersVapour ,round(LayersIce), round(LayersCH1),round(LayersCHEndZone), RadiusVapour, RadiusIce, RadiusCHEndZone,RadiusCH, DensityVapour, DensityIce, DensityCH);
    [ZoneMass, MassDiff] = MassDifference(Radii, Density);
    
    %Round the relative parameters
    LayersIce=round(LayersIce);
    LayersCH1=round(LayersCH1);
    LayersCHEndZone = round(LayersCHEndZone);
    RatioIce = round(RatioIce, 4);
    RatioCH1 = round(RatioCH1,4);
    RatioCHEndZone = round(RatioCHEndZone,4);
    
    %Set the unused parameters to nan;
    LayersCH = nan;
    RatioCH = nan;
    
    %Calculate the number of layers, and the mass difference.
    NumLayers = LayersCH1+LayersCHEndZone+LayersVapour+LayersIce;
    MaxDiff = max(abs(MassDiff(LayersVapour:end)));
    
    %If there are too few intermediate layers, the mass matching can be
    %poor giving high percentage mass difference. In this case, error so
    %that the other approach is tried.
    if LayersCH1<50
        error
    end
    
    %If manual mode is selected, plot.
    if Figures==1
        figure
        plot(MassDiff)
        ylabel('Percentage Mass Difference')
        yyaxis right;
        plot(ZoneMass, ':');
        ylabel('Zone Mass')
        title('Zoning Plot')
        xlabel('Zone')
        xline(BoundaryVapour+LayersIce);
        figure
        plot(Radii(BoundaryVapour:end-2), MassDiff(BoundaryVapour:end))
        ylabel('Percentage Mass Difference')
        yyaxis right;
        plot(Radii(BoundaryVapour:end-1), ZoneMass(BoundaryVapour:end), ':');
        ylabel('Zone Mass')
        title('Zoning Plot')
        xlabel('Zone')
        xline(Radii(BoundaryVapour+LayersIce));
        xline(Radii(BoundaryVapour+LayersIce+LayersCH1));
        figure
        plot(Radii(BoundaryVapour:end-2), MassDiff(BoundaryVapour:end))
        ylabel('Percentage Mass Difference')
        yyaxis right;
        plot(Radii(BoundaryVapour:end-1), ZoneMass(BoundaryVapour:end)./(4*pi()*Radii(BoundaryVapour+1:end).^2), ':');
        ylabel('Zone Mass')
        title('Zoning Plot')
        xlabel('Zone')
        xline(Radii(BoundaryVapour+LayersIce));
        xline(Radii(BoundaryVapour+LayersIce+LayersCH1));
    end
    
    %Method tells us which approach was used for the meshing (allows us to
    %identify errors). Regions tells us how many regions were used, for use
    %in other scripts.
    Method=1;
    Regions=4;
    
    %2) This method is the same as above, but at the split the thickness is
    %instead chosen by using the same relative thickness increase as is seen in
    %the ice layer. This should be a little bit more stable - it means we don't
    %reach a scenario where we ask more of the meshing for the CH than we asked
    %for the ice.
catch
    try
        
        ThicknessAtCHForIceUsingIce =4.5e+03*(RadiusIce-RadiusVapour)*(1/SizeMod) * ThicknessAtVapour *1/(DensityVapour/0.00065); %4.5e+03
        ThicknessAtCHForIce = min([ThicknessAtCHForIceUsingIce, ThicknessAtCHForIceUsingCH]);
        ThicknessAtCHForCH = ThicknessAtCHForIce*(DensityIce/DensityCH);
        
        %Solve for Ice Layer
        [RatioIce, LayersIce] = ThicknessSolver(RadiusVapour, RadiusIce, ThicknessAtVapour, ThicknessAtCHForIce);
        
        %Calculate thickness for CH end zone, and gradient from the ice region.
        %Use this to calculate radius that the end zone starts at
        ThicknessAtCHEndZone = (ThicknessAtCHForIce/ThicknessAtVapour)*ThicknessAtVacuum;
        Gradient = (ThicknessAtCHForIce - ThicknessAtVapour)/(RadiusIce - RadiusVapour);
        RadiusCHEndZone = RadiusCH - ((ThicknessAtCHEndZone - ThicknessAtVacuum)/Gradient);
        RadiusCHEndZone = round(RadiusCHEndZone,4);
        
        
        %Solve for the first CH zone
        [RatioCH1, LayersCH1] = ThicknessSolver(RadiusIce, RadiusCHEndZone, ThicknessAtCHForCH, ThicknessAtCHEndZone);
        %Solve for the end zone
        [RatioCHEndZone, LayersCHEndZone] = ThicknessSolver(RadiusCHEndZone, RadiusCH, ThicknessAtCHEndZone, ThicknessAtVacuum);
        
        %Combine all and plot
        [Radii, Density] = BoundarySolver4Layer(RatioVapour, round(RatioIce, 4), round(RatioCH1, 4),  round(RatioCHEndZone, 4),LayersVapour ,round(LayersIce), round(LayersCH1),round(LayersCHEndZone), RadiusVapour, RadiusIce, RadiusCHEndZone,RadiusCH, DensityVapour, DensityIce, DensityCH);
        [ZoneMass, MassDiff] = MassDifference(Radii, Density);
        
        %Round
        LayersIce=round(LayersIce);
        LayersCH1=round(LayersCH1);
        LayersCHEndZone = round(LayersCHEndZone);
        RatioIce = round(RatioIce, 4);
        RatioCH1 = round(RatioCH1,4);
        RatioCHEndZone = round(RatioCHEndZone,4);
        
        %Set the unused parameters to nan;
        LayersCH = nan;
        RatioCH = nan;
    
        %Calculate the number of layers, and the mass difference
        NumLayers = LayersCH1+LayersCHEndZone+LayersVapour+LayersIce;
        MaxDiff = max(abs(MassDiff(LayersVapour:end)));
        
        %If there are too few intermediate layers, the mass matching can be
        %poor giving high percentage mass difference. In this case, error so
        %that the other approach is tried.
        if LayersCH1<50
            error
        end
        
                      
        %If manual mode is selected, plot.
        if Figures==1
            figure
            plot(MassDiff)
            ylabel('Percentage Mass Difference')
            yyaxis right;
            plot(ZoneMass, ':');
            ylabel('Zone Mass')
            title('Zoning Plot')
            xlabel('Zone')
            xline(BoundaryVapour+LayersIce);
            
            figure
            plot(Radii(BoundaryVapour:end-2), MassDiff(BoundaryVapour:end))
            ylabel('Percentage Mass Difference')
            yyaxis right;
            plot(Radii(BoundaryVapour:end-1), ZoneMass(BoundaryVapour:end), ':');
            ylabel('Zone Mass')
            title('Zoning Plot')
            xlabel('Zone')
            xline(Radii(BoundaryVapour+LayersIce));
            xline(Radii(BoundaryVapour+LayersIce+LayersCH1));
            
            figure
            plot(Radii(BoundaryVapour:end-2), MassDiff(BoundaryVapour:end))
            ylabel('Percentage Mass Difference')
            yyaxis right;
            plot(Radii(BoundaryVapour:end-1), ZoneMass(BoundaryVapour:end)./(4*pi()*Radii(BoundaryVapour+1:end).^2), ':');
            ylabel('Zone Mass')
            title('Zoning Plot')
            xlabel('Zone')
            xline(Radii(BoundaryVapour+LayersIce));
            xline(Radii(BoundaryVapour+LayersIce+LayersCH1));
        end
        
        %Identify method and regions
        Method=2;
        Regions=4;
        
        %Finally, if both methods have failed, we return to the 3 layer solution.
        %This normally works if the thickness/percentage mass difference is set
        %low enough, but can produce large numbers of layers if the CH region is
        %large. However, if the previous two have errored this usually means that
        %the CH and ice regions are roughly equal in size, and so this gives a good
        %solution. I did try a reverse of the above solution for when the ice layer
        %is larger, but this did not work well. In reality, the density difference
        %means (I think) that this scenario is uncommon, and can be described well
        %using the 3 layer approach.
    catch
        
        ThicknessAtCHForIceUsingIce =4.5e+03*(RadiusIce-RadiusVapour)*(1/SizeMod) * ThicknessAtVapour *1/(DensityVapour/0.00065); %4.5e+03
        ThicknessAtCHForIceUsingCH = 15e+03 * (RadiusCH-RadiusIce) * ThicknessAtVacuum  / (DensityIce/DensityCH); %15e+03
        ThicknessAtCHForIce = min([ThicknessAtCHForIceUsingIce, ThicknessAtCHForIceUsingCH]);
        ThicknessAtCHForCH = ThicknessAtCHForIce*(DensityIce/DensityCH);
        
        %Solve for Ice Layer and CH Layer
        [RatioIce, LayersIce] = ThicknessSolver(RadiusVapour, RadiusIce, ThicknessAtVapour, ThicknessAtCHForIce);
        [RatioCH, LayersCH] = ThicknessSolver(RadiusIce, RadiusCH, ThicknessAtCHForCH, ThicknessAtVacuum);
        
        %Use these values to create Radii, Density etc. info.
        [Radii, Density] = BoundarySolver(RatioVapour, round(RatioIce, 4), round(RatioCH, 4), LayersVapour ,round(LayersIce), round(LayersCH), RadiusVapour, RadiusIce, RadiusCH, DensityVapour, DensityIce, DensityCH);
        [ZoneMass, MassDiff] = MassDifference(Radii, Density);
        
        %Round
        LayersIce = round(LayersIce);
        LayersCH = round(LayersCH);
        RatioIce = round(RatioIce, 4);
        RatioCH = round(RatioCH,4);
        
        %Set the unused parameters to nan;
        LayersCH1 = nan;
        RatioCH1 = nan;
        LayersCHEndZone = nan;
        RatioCHEndZone = nan;
        RadiusCHEndZone = nan;
    
        %Calculate number of layers, mass difference
        NumLayers = LayersCH+LayersVapour+LayersIce;
        MaxDiff = max(abs(MassDiff(LayersVapour:end)));
        
        %Plot if in manual mode
        if Figures==1
            figure
            plot(MassDiff)
            ylabel('Percentage Mass Difference')
            yyaxis right;
            plot(ZoneMass, ':');
            ylabel('Zone Mass')
            title('Zoning Plot')
            xlabel('Zone')
            xline(BoundaryVapour+LayersIce);
        end
        
        %Identify method and number of regions
        Method=3;
        Regions=3;
    end
end











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

%Given the radii and density, will calculate the zone masses and return
%the percentage mass difference
function [ZoneMass, MassDiff] = MassDifference(Radii, Density)
Vol = (4/3) * pi() * (Radii.^3);
ZoneMass = diff(Vol).*Density;
MassDiff = 100*(diff(ZoneMass)./ZoneMass(2:end));
end

%Optimisation function. Given the Ratios, Radii and Layers, will construct
%the mesh and calculate the mass difference using above functions. This
%solves for 3 layer systems
function [Radii, Density] = BoundarySolver(RatioVapour, RatioIce, RatioCH, LayersVapour,  LayersIce, LayersCH, RadiusVapour, RadiusIce, RadiusCH, DensityVapour, DensityIce, DensityCH)

BoundaryVapour = LayersVapour+1;
BoundaryIce = BoundaryVapour + LayersIce;
BoundaryCH = BoundaryIce + LayersCH;

VapourRadii = RatioIncrement(1, BoundaryVapour, 0, RadiusVapour, RatioVapour);
IceRadii = RatioIncrement(BoundaryVapour, BoundaryIce, RadiusVapour, RadiusIce, RatioIce);
CHRadii = RatioIncrement(BoundaryIce, BoundaryCH, RadiusIce, RadiusCH, RatioCH);

Density = [DensityVapour*ones(1, LayersVapour), DensityIce*ones(1, (LayersIce)), DensityCH*ones(1, LayersCH)];

Radii = [0, VapourRadii, IceRadii, CHRadii];


end

%Optimisation function. Given the Ratios, Radii and Layers, will construct
%the mesh and calculate the mass difference using above functions. This
%solves for 4 layers sytems where the split is in the CH.
function [Radii, Density] = BoundarySolver4Layer(RatioVapour, RatioIce, RatioCH1,RatioCH2, LayersVapour,  LayersIce, LayersCH1,LayersCH2, RadiusVapour, RadiusIce,RadiusSplit, RadiusCH, DensityVapour, DensityIce, DensityCH)

BoundaryVapour = LayersVapour+1;
BoundaryIce = BoundaryVapour + LayersIce;
BoundaryCH1 = BoundaryIce + LayersCH1;
BoundaryCH2 = BoundaryCH1 + LayersCH2;

VapourRadii = RatioIncrement(1, BoundaryVapour, 0, RadiusVapour, RatioVapour);
IceRadii = RatioIncrement(BoundaryVapour, BoundaryIce, RadiusVapour, RadiusIce, RatioIce);
CH1Radii = RatioIncrement(BoundaryIce, BoundaryCH1, RadiusIce, RadiusSplit, RatioCH1);
CH2Radii = RatioIncrement(BoundaryCH1, BoundaryCH2, RadiusSplit, RadiusCH, RatioCH2);

Density = [DensityVapour*ones(1, LayersVapour), DensityIce*ones(1, (LayersIce)), DensityCH*ones(1, LayersCH1), DensityCH*ones(1, LayersCH2)];

Radii = [0, VapourRadii, IceRadii, CH1Radii, CH2Radii];

end


%Optimisation function. Given the Ratios, Radii and Layers, will construct
%the mesh and calculate the mass difference using above functions. This
%solves for 4 layers sytems where the split is in the ice.
function [Radii, Density] = BoundarySolver4Layer2(RatioVapour, RatioIce1, RatioIce2,RatioCH, LayersVapour,  LayersIce1, LayersIce2,LayersCH, RadiusVapour, RadiusIce1,RadiusIce, RadiusCH, DensityVapour, DensityIce, DensityCH)

BoundaryVapour = LayersVapour+1;
BoundaryIce1 = BoundaryVapour + LayersIce1;
BoundaryIce = BoundaryIce1 + LayersIce2;
BoundaryCH = BoundaryIce + LayersCH;

VapourRadii = RatioIncrement(1, BoundaryVapour, 0, RadiusVapour, RatioVapour);
Ice1Radii = RatioIncrement(BoundaryVapour, BoundaryIce1, RadiusVapour, RadiusIce1, RatioIce1);
IceRadii = RatioIncrement(BoundaryIce1, BoundaryIce, RadiusIce1, RadiusIce, RatioIce2);
CHRadii = RatioIncrement(BoundaryIce, BoundaryCH, RadiusIce, RadiusCH, RatioCH);

Density = [DensityVapour*ones(1, LayersVapour), DensityIce*ones(1, (LayersIce1)), DensityIce*ones(1, LayersIce2), DensityCH*ones(1, LayersCH)];

Radii = [0, VapourRadii, Ice1Radii, IceRadii, CHRadii];

end


%Function to solve for optimal ratio and number of layers. Given the radii
%and desired thickness at either end of the region, this function will
%solve the two simultaneous equations to determine the optimal number of
%layers and ratio to use.
function [Ratio, Layers] = ThicknessSolver(LowerRadius, UpperRadius, LowerThickness, UpperThickness)
%Turn off numerical solver warning message
warning('off','symbolic:solve:FallbackToNumerical');

%Solve analytic formula to find ratio and layers for ice region .
syms a n
eqn1 = log(1 - (1-a)*(UpperRadius-LowerRadius)/LowerThickness)/log(a) == n;
eqn2 = a == (UpperThickness/LowerThickness)^(1/(n-1));
sol = solve([eqn1, eqn2], [a, n], 'Real', true);
aSol = sol.a;
nSol = sol.n;
Ratio=double(aSol);
Layers=double(nSol);

end

