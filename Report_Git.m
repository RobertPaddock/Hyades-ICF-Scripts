%Produces big report card with text and multiple plots. Called within Analysis script. Subplot puts
%figures too far together, so the position commands move the endmost plots.
figure

set(gcf, 'Position', get(0, 'Screensize'));
% set(gcf, 'Position',  [0, 0, 2160, 1080])
subplot(1,6, [1,2])
dim = [.1 .57 .25 .35];

strA = ['Vapour Radius = ', num2str(IceBoundary(1))];
strB = [', Ice Radius = ', num2str(CHBoundary(1))];
strC = [', CH Radius = ', num2str(InitialRadius), newline];
strD = ['Vapour Pressure = ', num2str(Density(1,1)), newline, newline];
if Type==1 
    strE = ['Pulse Powers = ', num2str((abs(PulsePowers)./10^12).', '%.2f TW, '), newline];
    strF = ['Estimated Switch On = ', num2str((abs(PulseOnTime)./10^-9).', '%.2f ns, '), newline];
else 
    strE = ' '; strF =' ';
end
strG=(['Power = ' ,num2str(max(TotalLaserPower)/10^12), ' TW, Energy = ' num2str(TotalLaserEnergy/1000), ' kJ, Required = ' num2str(RequiredLaserEnergy/1000), ' kJ.', newline,newline, 'Neutron energy = '  num2str(NeutronEnergy/1000) ' kJ, Gain = ' num2str(Gain) newline]);
% if max(OlsonCR)<20 || Type==0
%     strH=(['Olson CR of ',num2str(MaxOlsonCR), ' (Lindl CR of ',num2str(max(LindlCR(~isinf(LindlCR)))), ' )', newline] );
% else
%     [~, ShellMinIndex] = min(ShellRadius);
%     strH=(['Olson CR of ',num2str(MaxOlsonCR), newline , '(Lindl CR of ',num2str(max(LindlCR(1:MeaningfulCRRange))),  ' )', newline] );
% end
strH=(['Olson CR of ',num2str(MaxOlsonCR), ' (Lindl CR of ',num2str(MaxLindlCR), ' )', newline] );
strI=(['IFAR of ' ,num2str(IFARCraxton), ' (Should be under 30)', newline]);
strJ=(['Implosion Velocity of ' ,num2str(abs(min(ImplosionVelocity))),  ' km/s (Should be under 400)',newline, newline]);
strK=(['Bang Time = ' ,num2str((BangTime./10^-9).', '%.4f ns ')]);
strL=sprintf(', Neutrons produced = %e \n', max(Neutrons));
strM=sprintf('Parametric Limit of %.3e (Should be under 10^1^4) \n', ParametricLimit);
strN=(['Burn Up fraction = ' ,num2str(BurnUpFraction), '%', newline]);
strO = ['Shell K.E. = ', num2str(ShellKineticEnergyValue/1000, '%.3f'), ' kJ, Hydroefficiency = ', num2str(HydroEfficiencyValue, '%.2f'), '%, Total Efficiency = ', num2str(TotalEfficiencyValue, '%.2f'), '%' newline];
strQ=(['Adiabat (Goncharov defn.) = ' ,num2str(AdiabatGoncharovValue), newline]);
strP=strrep(file, '\', '\\');
if strP(1) == 'C'
strfile1 = strP(1:39);
strfile2 = strP(40:end);
strfile = [strfile1 newline strfile2];
elseif strP(1) == '\'
    strfile1 = strP(1:49);
strfile2 = strP(50:end);
strfile = [strfile1 newline strfile2];
else 
    strfile = strP;
end
str = [strA strB strC strD strE strF strG strH strI strJ strK strL strM strN strO strQ newline strfile];
annotation('textbox',dim,'String',str)

subplot(2,6,3)
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
pos = get(gca, 'Position');
pos(1) = pos(1) + 0.00;
set(gca, 'Position', pos)

subplot(2,6, 4)
plot(HSAverageTemp(MinStagnationIndex:MaxStagnationIndex), HSAverageOldRhoR(MinStagnationIndex:MaxStagnationIndex));
hold on
plot(SimT, SimOldRhoR)
xlabel('Average Hot Spot Ion Temp (eV)')
ylabel('Average Hot Spot OldRhoR (g/cm^2)')
legend ('Simulated data', 'Cheng analytic formula')
TimestepsAboveIgnition=sum(HSAverageOldRhoR>4*5.514*(HSAverageTemp/1000).^-2.5);
pos = get(gca, 'Position');
pos(1) = pos(1) + 0.01;
set(gca, 'Position', pos)

subplot(2,6,5)
title('Hotspot')
yyaxis left;
plot(Time, HSAverageTemp);
ylabel('Temperature, eV');
xlim([Time(MinStagnationIndex) Time(MaxStagnationIndex)])
hold on;
yyaxis right;
plot(Time, HSRhoR);
ylabel('RhoR, g/cm^2');
xlim([Time(MinStagnationIndex) Time(MaxStagnationIndex)])
xline(BangTime, 'k');
pos = get(gca, 'Position');
pos(1) = pos(1) + 0.01;
set(gca, 'Position', pos)

subplot(2,6,6)
title('Shell')
yyaxis left;
plot(Time, ShellAverageTemp);
ylabel('Temperature, eV');
xlim([Time(MinStagnationIndex) Time(MaxStagnationIndex)])
hold on;
yyaxis right;
plot(Time, ShellRhoR)
ylabel('RhoR, g/cm^2');
xlim([Time(MinStagnationIndex) Time(MaxStagnationIndex)])
xline(BangTime, 'k');
pos = get(gca, 'Position');
pos(1) = pos(1) + 0.03;
set(gca, 'Position', pos)

subplot(2,6, [7,8])
yyaxis left
plot(Time./10^-9, IFAR)
ylabel('IFAR');
yyaxis right
plot(Time./10^-9, abs(ImplosionVelocity))
ylabel('Velocity (km/s)');
xlabel('Time (ns)');
pos = get(gca, 'Position');
pos(1) = pos(1) - 0.02;
set(gca, 'Position', pos)

subplot(2,6, [9,10]);
plot(Time.*10^9, AdiabatGoncharov)
xlabel('Time (ns)')
ylabel('Adiabat (Goncharov definition)')
xlim([0 BangTime.*10^9+0.2])
xline(BangTime.*10^9)
xline(Time(AdiabatTimeIndex).*10^9)

subplot(2,6, [11,12]);
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
 pos = get(gca, 'Position');
 pos(1) = pos(1) + 0.02;
 set(gca, 'Position', pos)
