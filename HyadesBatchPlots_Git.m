%Plots data loaded using Hyades Batch Load

Mode = 3; %1 for thickness, 2 for laser (2.5 for 4 pulses), 3 for electron heating (3 for timing, 3.5 for power), 4 for cutoff, 5 for vapour pressure, 6 for second laser (6.5 for second laser cutoff)
Order = 1; %0 or 1, which is fixed variable and which is change
All =1; %1 for scatter of all data, 0 for plots following order
ShowAllPlots = 0;

%For the selected mode and order, saves the fixed and changing values for
%the valid files into the corresponding index (there will be more initial
%variables than files since some files error or are missing - using valid
%files connects the two up)
if Mode == 2
    Pulse2Time = PulseAll(2,:);
    Pulse3Time = PulseAll(3,:);
    Pulse2TimeValues=unique(Pulse2Time);
    Pulse3TimeValues=unique(Pulse3Time);
    if Order == 1
        FixedVariableValues = Pulse2TimeValues;
        FixedVariable = Pulse2Time;
         ChangeVariable = Pulse3Time;
         Xtitle = 'Pulse 3 Time';
         Legendtitle = 'Pulse 2 Time';
    else
        FixedVariableValues = Pulse3TimeValues;
        FixedVariable = Pulse3Time;
        ChangeVariable = Pulse2Time;
        Xtitle = 'Pulse 2 Time';
        Legendtitle = 'Pulse 3 Time';
    end
elseif Mode == 2.5
    Pulse2Time = PulseAll(2,:);
    Pulse3Time = PulseAll(3,:);
    Pulse4Time = PulseAll(4,:);
    Pulse2TimeValues=unique(Pulse2Time);
    Pulse3TimeValues=unique(Pulse3Time);
    Pulse4TimeValues=unique(Pulse4Time);
    if Order == 1
        FixedVariableValues = Pulse3TimeValues;
        FixedVariable = Pulse3Time;
         ChangeVariable = Pulse4Time;
         Xtitle = 'Pulse 4 Time';
         Legendtitle = 'Pulse 3 Time';
    else
        FixedVariableValues = Pulse4TimeValues;
        FixedVariable = Pulse4Time;
        ChangeVariable = Pulse3Time;
        Xtitle = 'Pulse 3 Time';
        Legendtitle = 'Pulse 4 Time';
    end
elseif Mode == 1
    RadiusIce = IceRadiusAll;
    RadiusCH = CHRadiusAll;
    RadiusVapour = VapourRadiusAll;
    RadiusIceValues = unique(RadiusIce);
    RadiusVapourValues = unique(RadiusVapour);
         if Order == 1
        FixedVariableValues = RadiusIceValues;
        FixedVariable = RadiusIce;
        ChangeVariable = RadiusVapour;
        Xtitle = 'Vapour/Ice Boundary';
        Legendtitle = 'Ice/CH Boundary';
    else
               FixedVariableValues = RadiusVapourValues;
        FixedVariable = RadiusVapour;
        ChangeVariable = RadiusIce;
        Legendtitle = 'Vapour/Ice Boundary';
        Xtitle = 'Ice/CH Boundary';
    end
elseif Mode==3
    ChangeVariable = HeatingTimingsAll(:,1).';
    FixedVariable = HeatingTimingsAll(:,1).';
    Xtitle = 'Electron Heating Applied';
    Legendtitle = 'Ice/CH Boundary';
    elseif Mode==3.5
    ChangeVariable = HeatingEnergyAll/1000;
    FixedVariable = HeatingEnergyAll/1000;
    Xtitle = 'Heating energy (kJ)';
    Legendtitle = 'Ice/CH Boundary';
 elseif Mode==4
    ChangeVariable = RequiredLaserEnergyAll;
    FixedVariable = RequiredLaserEnergyAll;
    Xtitle = 'Laser Energy';
    Legendtitle = 'Ice/CH Boundary';
 elseif Mode==5
    ChangeVariable = VapourPressureAll;
    FixedVariable = VapourPressureAll;
    Xtitle = 'Vapour Density';
    Legendtitle = 'Vapour Density';
elseif Mode==6
    ChangeVariable = SecondLaserOnTimeAll;
    FixedVariable = SecondLaserOnTimeAll;
    Xtitle = 'Second Laser Time';
    Legendtitle = 'Second Laser Time';
elseif Mode==6.5
    ChangeVariable = SecondLaserOffTimeAll;
    FixedVariable = SecondLaserOffTimeAll;
    Xtitle = 'Second Laser Off Time';
    Legendtitle = 'Second Laser Off Time';
end

map=[1 0 0; 0 0 1];
map2=[ 0 0 0; 0 0.5 0];
ColorCR = OlsonCRAll>16;
ImplosionVelocityLimit = 400;
ColorVelocity = (ImplosionVelocityAll<ImplosionVelocityLimit) .* (IFARAll<30);
 pointsize = 20;
if All==1 
    if ShowAllPlots ==1
    %Plot the data, with a datatip saying the file name
    fig1 = figure('DeleteFcn','datacursormode');
   title('Neutrons (consider end time!)')
    scatter(ChangeVariable,NeutronsAll,pointsize, ColorCR);
    xlabel(Xtitle);
    ylabel('Neutron Yield')
    colorbar()
    dcm_obj = datacursormode(fig1);
    set(dcm_obj,'UpdateFcn',{@myupdatefcn,ValidFiles,GainAll,OlsonCRAll, ChangeVariable, FixedVariable});
    colormap(map)
    caxis([0 1])

    fig2 = figure('DeleteFcn','datacursormode');
    title('Convergence Ratio')
    scatter(ChangeVariable,OlsonCRAll,pointsize, ColorCR);
    xlabel(Xtitle);
    ylabel('CR')
    colorbar()
    dcm_obj = datacursormode(fig2);
    set(dcm_obj,'UpdateFcn',{@myupdatefcn,ValidFiles,GainAll,OlsonCRAll, ChangeVariable, FixedVariable});
    colormap(map)
    caxis([0 1])

    fig3 = figure('DeleteFcn','datacursormode');
    title('IFAR')
    scatter(ChangeVariable,IFARAll,pointsize, ColorCR);
    xlabel(Xtitle);
    ylabel('IFAR')
    colorbar()
    dcm_obj = datacursormode(fig3);
    set(dcm_obj,'UpdateFcn',{@myupdatefcn,ValidFiles,GainAll,OlsonCRAll, ChangeVariable, FixedVariable});
    colormap(map)
    caxis([0 1])

    fig4 = figure('DeleteFcn','datacursormode');
    title('Implosion Velocity')
    scatter(ChangeVariable,ImplosionVelocityAll,pointsize, ColorCR);
    xlabel(Xtitle);
    ylabel('Implosion Velocity')
    colorbar()
    dcm_obj = datacursormode(fig4);
    set(dcm_obj,'UpdateFcn',{@myupdatefcn,ValidFiles,GainAll,OlsonCRAll, ChangeVariable, FixedVariable});
    colormap(map)
    caxis([0 1])

    fig5 = figure('DeleteFcn','datacursormode');
    title('Gain')
    scatter(ChangeVariable,GainAll,pointsize, ColorCR);
    xlabel(Xtitle);
    ylabel('Gain')
    colorbar()
    dcm_obj = datacursormode(fig5);
    set(dcm_obj,'UpdateFcn',{@myupdatefcn,ValidFiles,GainAll,OlsonCRAll, ChangeVariable, FixedVariable})
    colormap(map)
    caxis([0 1])

    fig6 = figure('DeleteFcn','datacursormode');
    title('Energy')
    scatter(ChangeVariable,EnergyAll,pointsize, ColorCR);
    xlabel(Xtitle);
    ylabel('Input Energy')
    colorbar()
    dcm_obj = datacursormode(fig6);
    set(dcm_obj,'UpdateFcn',{@myupdatefcn,ValidFiles,GainAll,OlsonCRAll, ChangeVariable, FixedVariable});
    colormap(map)
    caxis([0 1])

    fig7 = figure('DeleteFcn','datacursormode');
    title('Time')
    scatter(ValidFiles,TimeAll,pointsize, ColorCR);
    xlabel(Xtitle);
    ylabel('End Time')
    colorbar()
    dcm_obj = datacursormode(fig7);
    set(dcm_obj,'UpdateFcn',{@myupdatefcn,ValidFiles,GainAll,OlsonCRAll, ChangeVariable, FixedVariable});
    colormap(map)
    caxis([0 1])

    fig8 = figure('DeleteFcn','datacursormode');
    title('Gain')
    scatter(OlsonCRAll,GainAll,pointsize, ColorCR);
    xlabel('CR')
    ylabel('Gain')
    colorbar()
    dcm_obj = datacursormode(fig8);
    set(dcm_obj,'UpdateFcn',{@myupdatefcn,ValidFiles,GainAll,OlsonCRAll, ChangeVariable, FixedVariable})
    colormap(map)
    caxis([0 1])

    fig9 = figure('DeleteFcn','datacursormode');
    title('Burn Up Fraction')
    scatter(ChangeVariable,BurnUpAll,pointsize, ColorCR);
    xlabel(Xtitle);
    ylabel('Fraction')
    colorbar()
    dcm_obj = datacursormode(fig9);
    set(dcm_obj,'UpdateFcn',{@myupdatefcn,ValidFiles,GainAll,OlsonCRAll, ChangeVariable, FixedVariable})
    colormap(map)
    caxis([0 1])
    
     fig10 = figure('DeleteFcn','datacursormode');
    title('Gain vs implosion velocity')
    scatter(ImplosionVelocityAll,GainAll,pointsize, ColorCR);
    xlabel('Implosion Velocity');
    ylabel('Gain')
    colorbar()
    dcm_obj = datacursormode(fig10);
    set(dcm_obj,'UpdateFcn',{@myupdatefcn,ValidFiles,GainAll,OlsonCRAll, ChangeVariable, FixedVariable})
    colormap(map)
    caxis([0 1])

    fig11 = figure('DeleteFcn','datacursormode');
    title('Gain')
    scatter(OlsonCRAll,GainAll,pointsize, ColorVelocity);
    xlabel('CR')
    ylabel('Gain')
    colorbar()
    dcm_obj = datacursormode(fig11);
    set(dcm_obj,'UpdateFcn',{@myupdatefcn,ValidFiles,GainAll,OlsonCRAll, ChangeVariable, FixedVariable})
    colormap(map)
    caxis([0 1])
      
    end
    %Multiplot showing important parameters
    figure
    set(gcf, 'Position',  [0, 0, 2160, 1080])
    subplot(2,3, 1)
    title('Gain')
    scatter(ChangeVariable,GainAll,pointsize, ColorCR);
    xlabel(Xtitle);
    ylabel('Gain')
    colormap(gca, map)
    caxis([0 1])

    subplot(2,3, 2)
    title('Implosion Velocity')
    scatter(ChangeVariable,ImplosionVelocityAll,pointsize, ColorCR);
    xlabel(Xtitle);
    ylabel('Implosion Velocity')
    colormap(gca,map)
    caxis([0 1])

    subplot(2,3, 3)
    title('Gain')
    scatter(OlsonCRAll,GainAll,pointsize, ColorVelocity);
    xlabel('CR')
    ylabel('Gain')
    colormap(gca,map2)
    caxis([0 1])

    subplot(2,3, 4)
    title('Convergence Ratio')
    scatter(ChangeVariable,OlsonCRAll,pointsize, ColorCR);
    xlabel(Xtitle);
    ylabel('CR')
    colormap(gca,map)
    caxis([0 1])

    subplot(2,3, 5)
    title('IFAR')
    scatter(ChangeVariable,IFARAll,pointsize, ColorCR);
    xlabel(Xtitle);
    ylabel('IFAR')
    colormap(gca,map)
    caxis([0 1])

    subplot(2,3, 6)
    title('Gain')
    scatter(ImplosionVelocityAll,GainAll,pointsize, ColorCR);
    xlabel('Implosion Velocity');
    ylabel('Gain')
    colormap(gca,map)
    caxis([0 1])

    dcm_obj = datacursormode(gcf);
    set(dcm_obj,'UpdateFcn',{@myupdatefcn,ValidFiles,GainAll,OlsonCRAll, ChangeVariable, FixedVariable});

    %Multiplot showing temps and densities
    figure
    set(gcf, 'Position',  [0, 0, 2160, 1080])
    subplot(2,3, 1)
    title('HS RhoR')
    scatter(ChangeVariable,MaxHSRhoRAll,pointsize, ColorCR);
    xlabel(Xtitle);
    ylabel('HSRhoR')
    colormap(map)
    caxis([0 1])
    subplot(2,3, 2)
    title('Shell RhoR')
    scatter(ChangeVariable,MaxShellRhoRAll,pointsize, ColorCR);
    xlabel(Xtitle);
    ylabel('ShellRhoR')
    colormap(map)
    caxis([0 1])
    subplot(2,3, 4)
    title('HS Temp')
    scatter(ChangeVariable,MaxHSTempAll,pointsize, ColorCR);
    xlabel(Xtitle);
    ylabel('HSTemp')
    colormap(map)
    caxis([0 1])
    subplot(2,3, 5)
    title('Shell Temp')
    scatter(ChangeVariable,MaxShellTempAll,pointsize, ColorCR);
    xlabel(Xtitle);
    ylabel('ShellTemp')
    colormap(map)
    caxis([0 1])
    subplot(2,3, 3)
    title('Burn Up Fraction')
    scatter(ChangeVariable,BurnUpAll,pointsize, ColorCR);
    xlabel(Xtitle);
    ylabel('Fraction')
    caxis([0 1])
    subplot(2,3, 6)
    title('Time')
    scatter(ValidFiles,TimeAll,pointsize, ColorCR);
    xlabel(Xtitle);
    ylabel('End Time')
    colormap(gca,map)
    caxis([0 1])

    dcm_obj = datacursormode(gcf);
    set(dcm_obj,'UpdateFcn',{@myupdatefcn,ValidFiles,GainAll,OlsonCRAll, ChangeVariable, FixedVariable});
else
    
    
    
    
    
    
%Sets figures up for use with custom datatips
fig1 = figure('DeleteFcn','datacursormode');
set(gcf, 'Position',  [0, 0, 2160, 1080]);
fig2 = figure('DeleteFcn','datacursormode');
set(gcf, 'Position',  [0, 0, 2160, 1080]);
if ShowAllPlots==1
fig7 = figure('DeleteFcn','datacursormode');
fig8 = figure('DeleteFcn','datacursormode');
fig3 = figure('DeleteFcn','datacursormode');
fig4 = figure('DeleteFcn','datacursormode');
fig5 = figure('DeleteFcn','datacursormode');
fig6 = figure('DeleteFcn','datacursormode');
end


%Pulse3TimeValues = [7:0.2:9]*10^-9;

    
%Chooses most distinguishable plot colours.
c = distinguishable_colors(length(FixedVariableValues));
i=1;

%Criteria identifies the relevant files that meet the criteria. We then
%create new variables corresponding only to the file that meets the
%criteria, and then plot these. By using the fixed variable values as
%criteria, we get a series of plots for each fixed variable value.
for Index = FixedVariableValues
    Criteria = FixedVariable==Index;

    NeutronPlot = NeutronsAll(Criteria);
    CRPlot = OlsonCRAll(Criteria);
    IFARPlot = IFARAll(Criteria);
    VelocityPlot = ImplosionVelocityAll(Criteria);
    GainPlot=GainAll(Criteria);
    PlotVariable = ChangeVariable(Criteria);
    MaxHSTempPlot = MaxHSTempAll(Criteria);
    MaxHSRhoRPlot = MaxHSRhoRAll(Criteria);
    MaxShellTempPlot = MaxShellTempAll(Criteria);
    MaxShellRhoRPlot = MaxShellRhoRAll(Criteria);
    BurnUpPlot = BurnUpAll(Criteria);
      
    if ShowAllPlots==1
    %Plot the figures
    figure(7)
    hold on
    plot(PlotVariable, NeutronPlot, 'color', c(i,:))
    title('Neutron Yield')
           
    figure(8)
    hold on
    plot(PlotVariable, CRPlot, 'color', c(i,:))
    title('Convergence Ratio')
    
    figure(3)
    hold on
    plot(PlotVariable, IFARPlot, 'color', c(i,:))
    title('IFAR')
    
     figure(4)
    hold on
    plot(PlotVariable, VelocityPlot, 'color', c(i,:))
    title('Peak Implosion Velocity')
    %dcm_obj = datacursormode(fig4);
    %set(dcm_obj,'UpdateFcn',{@myupdatefcn,ValidFiles,GainAll, OlsonCRAll, ChangeVariable, FixedVariable});

    figure(5)
    hold on
    plot(PlotVariable, GainPlot, 'color', c(i,:))
    title('Gain')
    %dcm_obj = datacursormode(fig5);
    %set(dcm_obj,'UpdateFcn',{@myupdatefcn,ValidFiles,GainAll, OlsonCRAll, ChangeVariable, FixedVariable});
    
    figure(6)
    hold on
    plot(CRPlot, GainPlot, 'color', c(i,:))
    title('Gain')
    %dcm_obj = datacursormode(fig5);
    %set(dcm_obj,'UpdateFcn',{@myupdatefcn,ValidFiles,GainAll, OlsonCRAll, ChangeVariable, FixedVariable});
    
    end
    
    figure(1)
    subplot(2,3, 1)
    hold on
    plot(PlotVariable, GainPlot, 'color', c(i,:))
    title('Gain')
    xlabel(Xtitle);
    ylabel('Gain');

    subplot(2,3, 2)
    hold on
    plot(PlotVariable, VelocityPlot, 'color', c(i,:))
    title('Peak Implosion Velocity')
    xlabel(Xtitle);
    ylabel('Implosion Velocity')

    subplot(2,3, 3)
    title('Gain vs CR')
    hold on
    plot(CRPlot, GainPlot, 'color', c(i,:))
    xlabel('CR')
    ylabel('Gain')

    subplot(2,3, 4)
    title('Convergence Ratio')
    hold on
    plot(PlotVariable, CRPlot, 'color', c(i,:))
    xlabel(Xtitle);
    ylabel('CR')

    subplot(2,3, 5)
    title('IFAR')
    hold on
    plot(PlotVariable, IFARPlot, 'color', c(i,:))
    xlabel(Xtitle);
    ylabel('IFAR')

dcm_obj = datacursormode(fig1);
set(dcm_obj,'UpdateFcn',{{@myupdatefcn,ValidFiles,GainAll,OlsonCRAll, ChangeVariable, FixedVariable}});
    
figure(2)
subplot(2,3, 1)
title('HS RhoR')
hold on
plot(PlotVariable, MaxHSRhoRPlot, 'color', c(i,:))
xlabel(Xtitle);
ylabel('HSRhoR')
colormap(map)
subplot(2,3, 2)
title('Shell RhoR')
hold on
plot(PlotVariable, MaxShellRhoRPlot, 'color', c(i,:))
xlabel(Xtitle);
ylabel('ShellRhoR')
colormap(map)
subplot(2,3, 4)
title('HS Temp')
hold on
plot(PlotVariable, MaxHSTempPlot, 'color', c(i,:))
xlabel(Xtitle);
ylabel('HSTemp')
colormap(map)
subplot(2,3, 5)
title('Shell Temp')
hold on
plot(PlotVariable, MaxShellTempPlot, 'color', c(i,:))
xlabel(Xtitle);
ylabel('ShellTemp')
colormap(map)
subplot(2,3, 3)
title('Burn Up Fraction')
hold on
plot(PlotVariable, BurnUpPlot, 'color', c(i,:))
xlabel(Xtitle);
ylabel('Fraction')
dcm_obj = datacursormode(fig2);
set(dcm_obj,'UpdateFcn',{@myupdatefcn,ValidFiles,GainAll,OlsonCRAll, ChangeVariable, FixedVariable});
    i=i+1;
end

if ShowAllPlots==1
    
%Plot the data, with a datatip saying the file name
figure(7);
pointsize = 20;
title('Neutrons (consider end time!)')
plt1=scatter(ChangeVariable,NeutronsAll,pointsize, ColorCR);
ylabel('Neutron Yield')
dcm_obj = datacursormode(fig7);
set(dcm_obj,'UpdateFcn',{@myupdatefcn,ValidFiles,GainAll, OlsonCRAll, ChangeVariable, FixedVariable});
colormap(map)
leg7=legend(num2str(FixedVariableValues.'), 'Location','EastOutside');
title(leg7, Legendtitle);
xlabel(Xtitle);

figure(8);
title('Convergence Ratio')
scatter(ChangeVariable,OlsonCRAll,pointsize, ColorCR);
ylabel('CR')
dcm_obj = datacursormode(fig8);
set(dcm_obj,'UpdateFcn',{@myupdatefcn,ValidFiles,GainAll, OlsonCRAll, ChangeVariable, FixedVariable});
colormap(map)
leg8=legend(num2str(FixedVariableValues.'), 'Location','EastOutside');
title(leg8, Legendtitle);
xlabel(Xtitle);

figure(3);
title('IFAR')
scatter(ChangeVariable,IFARAll,pointsize, ColorCR);
ylabel('IFAR')
dcm_obj = datacursormode(fig3);
set(dcm_obj,'UpdateFcn',{@myupdatefcn,ValidFiles,GainAll, OlsonCRAll, ChangeVariable, FixedVariable});
colormap(map)
leg3=legend(num2str(FixedVariableValues.'), 'Location','EastOutside');
title(leg3, Legendtitle);
xlabel(Xtitle);

figure(4);
title('Implosion Velocity')
scatter(ChangeVariable,ImplosionVelocityAll,pointsize, ColorCR);
ylabel('Implosion Velocity')
dcm_obj = datacursormode(fig4);
set(dcm_obj,'UpdateFcn',{@myupdatefcn,ValidFiles,GainAll, OlsonCRAll, ChangeVariable, FixedVariable});
colormap(map)
leg4=legend(num2str(FixedVariableValues.'), 'Location','EastOutside');
title(leg4, Legendtitle);
xlabel(Xtitle);

figure(5);
title('Gain')
scatter(ChangeVariable,GainAll,pointsize, ColorCR);
ylabel('Gain')
dcm_obj = datacursormode(fig5);
set(dcm_obj,'UpdateFcn',{@myupdatefcn,ValidFiles,GainAll, OlsonCRAll, ChangeVariable, FixedVariable})
colormap(map)
leg5=legend(num2str(FixedVariableValues.'), 'Location','EastOutside');
title(leg5, Legendtitle);
xlabel(Xtitle);

figure(6);
title('Gain vs CR')
scatter(OlsonCRAll,GainAll,pointsize, ColorCR);
ylabel('Gain')
dcm_obj = datacursormode(fig6);
set(dcm_obj,'UpdateFcn',{@myupdatefcn,ValidFiles,GainAll, OlsonCRAll, ChangeVariable, FixedVariable})
colormap(map)
leg6=legend(num2str(FixedVariableValues.'), 'Location','EastOutside');
title(leg6, Legendtitle);
xlabel('CR');

end

figure(1)
leg1=legend(num2str(FixedVariableValues.'), 'Position',[0.765238144591443 0.178672556061061 0.0665280646692938 0.219960271897832]);
title(leg1, Legendtitle);

figure(2)
leg2=legend(num2str(FixedVariableValues.'), 'Position',[0.765238144591443 0.178672556061061 0.0665280646692938 0.219960271897832]);
title(leg2, Legendtitle);
end

function txt = myupdatefcn(~,event_obj,ValidFiles,GainAll, OlsonCRAll, ChangeVariable, FixedVariable)
% Customizes text of data tips
pos = get(event_obj,'Position');
I = get(event_obj, 'DataIndex');
txt = {['X: ',num2str(pos(1))],...
       ['Y: ',num2str(pos(2))],...
       ['Gain: ',num2str(GainAll(I))],...
       ['CR: ',num2str(OlsonCRAll(I))],...
       ['Change Variable: ',num2str(ChangeVariable(I))],...
       ['Fixed Variable: ',num2str(FixedVariable(I))],...
       ['File',num2str(ValidFiles(I))]};
end