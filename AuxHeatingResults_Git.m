ConversionEfficiency = 0.18*0.52;
ConversionAdjustment = 1/ConversionEfficiency;

%% Open all my heating files
%Each results.mat file has been copied to another folder, and renamed to
%refer to what the simulations used to produce it were investigating. These
%results files are then fed into a function (defined in this script) which
%extracts/calculates the key data (so this can be plotted later in the
%script).
file = ('\\aldaq1.physics.ox.ac.uk\Archer\Robert\Desktop Hyades Files\Dated\210617\210618ElectronHeating\Power0pt5.mat')
[HeatingEnergyPt5, LaserEnergyPt5, GainPt5, YieldPt5, TotalInputEnergyPt5, GainWithHeatingPt5, BurnUpPt5] = MassLoad(file, ConversionAdjustment);

file = ('\\aldaq1.physics.ox.ac.uk\Archer\Robert\Desktop Hyades Files\Dated\210617\210618ElectronHeating\Power0pt65.mat')
[HeatingEnergyPt65, LaserEnergyPt65, GainPt65, YieldPt65, TotalInputEnergyPt65, GainWithHeatingPt65, BurnUpPt65] = MassLoad(file, ConversionAdjustment);

file = ('\\aldaq1.physics.ox.ac.uk\Archer\Robert\Desktop Hyades Files\Dated\210617\210618ElectronHeating\Power0pt25.mat')
[HeatingEnergyPt25, LaserEnergyPt25, GainPt25, YieldPt25, TotalInputEnergyPt25, GainWithHeatingPt25, BurnUpPt25] = MassLoad(file, ConversionAdjustment);

file = ('\\aldaq1.physics.ox.ac.uk\Archer\Robert\Desktop Hyades Files\Dated\210617\210618ElectronHeating\Power0pt35.mat')
[HeatingEnergyPt35, LaserEnergyPt35, GainPt35, YieldPt35, TotalInputEnergyPt35, GainWithHeatingPt35, BurnUpPt35] = MassLoad(file, ConversionAdjustment);


%% Open the new heating files
file = '\\aldaq1.physics.ox.ac.uk\Archer\Robert\Tat summer project\TwoColour_0.65_Power_Input\Results.mat';
[HeatingEnergyTwoColourPt65, LaserEnergyTwoColourPt65, GainTwoColourPt65, YieldTwoColourPt65, TotalInputEnergyTwoColourPt65, GainWithHeatingTwoColourPt65, BurnUpTwoColourPt65] = MassLoad(file, ConversionAdjustment);
% figure; scatter(TotalInputEnergyTwoColourPt65, GainWithHeatingTwoColourPt65);title('Two Colour 0.65 Size'); xlabel('Total Energy'); ylabel('Gain')
% figure; scatter(TotalInputEnergyTwoColourPt65, BurnUpTwoColourPt65);title('Two Colour 0.65 Size'); xlabel('Total Energy'); ylabel('Burn up fraction')

file = '\\aldaq1.physics.ox.ac.uk\Archer\Robert\Tat summer project\TwoColour_0.5_Power\Results.mat';
[HeatingEnergyTwoColourPt5, LaserEnergyTwoColourPt5, GainTwoColourPt5, YieldTwoColourPt5, TotalInputEnergyTwoColourPt5, GainWithHeatingTwoColourPt5, BurnUpTwoColourPt5] = MassLoad(file, ConversionAdjustment);
% figure; scatter(TotalInputEnergyTwoColourPt5, GainWithHeatingTwoColourPt5); title('Two Colour 0.5 Size'); xlabel('Total Energy'); ylabel('Gain')
% figure; scatter(TotalInputEnergyTwoColourPt5, BurnUpTwoColourPt5); title('Two Colour 0.5 Size'); xlabel('Total Energy'); ylabel('Burn up fraction')
% % figure; scatter(HeatingEnergyTwoColourPt5, GainWithHeatingTwoColourPt5); title('Two Colour 0.5 Size'); xlabel('Heating Energy'); ylabel('Gain')

file = '\\aldaq1.physics.ox.ac.uk\Archer\Robert\Tat summer project\ArF_0.5_Power_60kJ_Input\Results.mat';
[HeatingEnergyArFPt5, LaserEnergyArFPt5, GainArFPt5, YieldArFPt5, TotalInputEnergyArFPt5, GainWithHeatingArFPt5, BurnUpArFPt5] = MassLoad(file, ConversionAdjustment);
% figure; scatter(TotalInputEnergyArFPt5, GainWithHeatingArFPt5); title('ArF 0.5 Size'); xlabel('Total Energy'); ylabel('Gain')
% figure; scatter(TotalInputEnergyArFPt5, BurnUpArFPt5); title('ArF 0.5 Size'); xlabel('Total Energy'); ylabel('Burn up fraction')


%% Open the RhoR heating files
file = '\\aldaq1.physics.ox.ac.uk\Archer\Robert\Tat summer project\Third Harmonic Areal Density\HeatingPower_Input\Results.mat';
[HeatingEnergyPt5RhoR, LaserEnergyPt5RhoR, GainPt5RhoR, YieldPt5RhoR, TotalInputEnergyPt5RhoR, GainWithHeatingPt5RhoR, BurnUpPt5RhoR] = MassLoad(file, ConversionAdjustment);

% figure; scatter(TotalInputEnergyPt5, GainWithHeatingPt5)
% hold on
% scatter(TotalInputEnergyPt5RhoR, GainWithHeatingPt5RhoR)
% hold off
% 
% figure; scatter(HeatingEnergyPt5, GainWithHeatingPt5)
% hold on
% scatter(HeatingEnergyPt5RhoR, GainWithHeatingPt5RhoR)
% hold off


file = '\\aldaq1.physics.ox.ac.uk\Archer\Robert\Tat summer project\Test\ArFOptimisedArealDensity\HeatPower_Output_File11\Results.mat';
[HeatingEnergyArFPt5RhoR, LaserEnergyArFPt5RhoR, GainArFPt5RhoR, YieldArFPt5RhoR, TotalInputEnergyArFPt5RhoR, GainWithHeatingArFPt5RhoR, BurnUpArFPt5RhoR] = MassLoad(file, ConversionAdjustment);
% figure; scatter(TotalInputEnergyArFPt5RhoR, GainWithHeatingArFPt5RhoR)
% hold on
% scatter(TotalInputEnergyArFPt5, GainWithHeatingArFPt5)
% hold off
% title('ArF optimised RhoR')
% 
% figure; scatter(HeatingEnergyArFPt5RhoR, GainWithHeatingArFPt5RhoR)
% hold on
% scatter(HeatingEnergyArFPt5, GainWithHeatingArFPt5)
% hold off
% title('ArF optimised RhoR')







%% Plot the new heating file data
cm = parula(20);

figure
subplot(2,2,1)
scatter(HeatingEnergyTwoColourPt5/1000, GainWithHeatingTwoColourPt5,  'ok', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
hold on
scatter(HeatingEnergyTwoColourPt65/1000, GainWithHeatingTwoColourPt65,'s',  'MarkerFaceColor', cm(11,:)', 'MarkerEdgeColor', 'k' )
scatter(HeatingEnergyArFPt5/1000, GainWithHeatingArFPt5, '^', 'MarkerFaceColor', cm(2,:)', 'MarkerEdgeColor', 'k')
hold off
xlabel({'Heating energy (kJ)'})
ylabel('Gain')
ylim([6 20])
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
fig=gcf;
fig.Units               = 'points';
fig.Position(3)         = 400;
fig.Position(4)         = 200;
set(fig.Children, ...
    'FontName',     'Times', ...
    'FontSize',     9);

subplot(2,2,2)
scatter(HeatingEnergyTwoColourPt5/1000, BurnUpTwoColourPt5, 'ok', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
hold on
scatter(HeatingEnergyTwoColourPt65/1000, BurnUpTwoColourPt65, 's',  'MarkerFaceColor', cm(11,:)', 'MarkerEdgeColor', 'k' )
scatter(HeatingEnergyArFPt5/1000, BurnUpArFPt5, '^', 'MarkerFaceColor', cm(2,:)', 'MarkerEdgeColor', 'k')
hold off
xlabel({'Heating energy (kJ)'})
ylabel('Burn up fraction')
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
fig=gcf;
fig.Units               = 'points';
fig.Position(3)         = 400;
fig.Position(4)         = 200;
set(fig.Children, ...
    'FontName',     'Times', ...
    'FontSize',     9);

subplot(2,2,[3 4])
scatter(TotalInputEnergyTwoColourPt5/1000, GainWithHeatingTwoColourPt5,  'ok', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
hold on
scatter(TotalInputEnergyTwoColourPt65/1000, GainWithHeatingTwoColourPt65, 's',  'MarkerFaceColor', cm(11,:)', 'MarkerEdgeColor', 'k')
scatter(TotalInputEnergyArFPt5/1000, GainWithHeatingArFPt5, '^', 'MarkerFaceColor', cm(2,:)', 'MarkerEdgeColor', 'k')
hold off
xlabel({'Total laser energy (kJ)'})
ylabel('Gain')
ylim([6 20])
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
set(fig.Children, ...
    'FontName',     'Times', ...
    'FontSize',     9);
fig=gcf;
fig.Units               = 'points';
fig.Position(3)         = 400;
fig.Position(4)         = 300;

annotation('textbox',[0.04 0.65 0.01 0.3],'String',{'(a)'},'FitBoxToText','on', 'EdgeColor', 'none');
annotation('textbox',[0.48 0.65 0.01 0.3],'String',{'(b)'},'FitBoxToText','on', 'EdgeColor', 'none');
annotation('textbox',[0.04 0.13 0.01 0.3],'String',{'(c)'},'FitBoxToText','on', 'EdgeColor', 'none');

annotation('textbox',[0.18 0.7025 0.06 0.07],'String',{'E'},'FitBoxToText','on', 'EdgeColor', 'none', 'Color', 'k');
annotation('textbox',[0.42,0.855,0.056,0.07],'String',{'F'},'FitBoxToText','on', 'EdgeColor', 'none', 'Color',	cm(11,:));
annotation('textbox',[0.29 0.785 0.06 0.07],'String',{'G'},'FitBoxToText','on', 'EdgeColor', 'none', 'Color',cm(2,:));

annotation('textbox',[0.64 0.72 0.06 0.07],'String',{'E'},'FitBoxToText','on', 'EdgeColor', 'none', 'Color', 'k');
annotation('textbox',[0.74,0.615,0.056,0.07],'String',{'F'},'FitBoxToText','on', 'EdgeColor', 'none', 'Color',	cm(11,:));
annotation('textbox',[0.82 0.855 0.06 0.07],'String',{'G'},'FitBoxToText','on', 'EdgeColor', 'none', 'Color',cm(2,:));

annotation('textbox',[0.21 0.21 0.06 0.07],'String',{'E'},'FitBoxToText','on', 'EdgeColor', 'none', 'Color', 'k');
annotation('textbox',[0.81,0.3875,0.056,0.07],'String',{'F'},'FitBoxToText','on', 'EdgeColor', 'none', 'Color',	cm(11,:));
annotation('textbox',[0.56 0.30 0.06 0.07],'String',{'G'},'FitBoxToText','on', 'EdgeColor', 'none', 'Color',cm(2,:));





%% Plot the RhoR data
cm = parula(20);

figure
subplot(1,2,1)
scatter(TotalInputEnergyPt5/1000, GainWithHeatingPt5,  'ok', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
hold on
scatter(TotalInputEnergyPt5RhoR/1000, GainWithHeatingPt5RhoR,'s',  'MarkerFaceColor', cm(11,:)', 'MarkerEdgeColor', 'k' )
hold off
xlabel({'Total input energy (kJ)'})
ylabel('Gain')
title('Third-harmonic 0.5 size')
xlim([700 1500])
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
fig=gcf;
fig.Units               = 'points';
fig.Position(3)         = 400;
fig.Position(4)         = 200;
set(fig.Children, ...
    'FontName',     'Times', ...
    'FontSize',     9);

subplot(1,2,2)
scatter(TotalInputEnergyArFPt5/1000, GainWithHeatingArFPt5, 'ok', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
hold on
scatter(TotalInputEnergyArFPt5RhoR/1000, GainWithHeatingArFPt5RhoR, 's',  'MarkerFaceColor', cm(11,:)', 'MarkerEdgeColor', 'k' )
hold off
xlabel({'Total input energy (kJ)'})
ylabel('Gain')
xlim([1900 2800])
title('ArF 0.5 size')
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
fig=gcf;
fig.Units               = 'points';
fig.Position(3)         = 400;
fig.Position(4)         = 200;
set(fig.Children, ...
    'FontName',     'Times', ...
    'FontSize',     9);

fig=gcf;
fig.Units               = 'points';
fig.Position(3)         = 400;
fig.Position(4)         = 200;

annotation('textbox',[0.04 0.65 0.01 0.3],'String',{'(a)'},'FitBoxToText','on', 'EdgeColor', 'none');
annotation('textbox',[0.48 0.65 0.01 0.3],'String',{'(b)'},'FitBoxToText','on', 'EdgeColor', 'none');

annotation('textbox',[0.164635320822472,0.493941951137357,0.144465287060496,0.157303367214703],'String',{'Optimised',  'Gain'},'FitBoxToText','on', 'EdgeColor', 'none', 'Color', 'k', 'HorizontalAlignment', 'center', 'FontName',     'Times');
annotation('textbox',[0.276392356469752,0.251923224545597,0.144465287060496,0.157303367214703],'String',{'Optimised',  '\rhoR'},'FitBoxToText','on', 'EdgeColor', 'none', 'Color',	cm(11,:), 'HorizontalAlignment', 'center', 'FontName',     'Times');

annotation('textbox',[0.6572,0.3037,0.1444,0.1575],'String',{'Optimised',  'Gain'},'FitBoxToText','on', 'EdgeColor', 'none', 'Color', 'k', 'HorizontalAlignment', 'center', 'FontName',     'Times');
annotation('textbox',[0.7658,0.7214,0.1445,0.1573],'String',{'Optimised',  '\rhoR'},'FitBoxToText','on', 'EdgeColor', 'none', 'Color',	cm(11,:), 'HorizontalAlignment', 'center', 'FontName',     'Times');



%% Plot overall

figure
scatter(HeatingEnergyPt25, YieldPt25./YieldPt25(1))
hold on
scatter(HeatingEnergyPt35, YieldPt35./YieldPt35(1))
scatter(HeatingEnergyPt5, YieldPt5./YieldPt5(1))
scatter(HeatingEnergyPt65, YieldPt65./YieldPt65(1))
scatter(HeatingEnergyTwoColourPt5, YieldTwoColourPt5./YieldTwoColourPt5(1))
scatter(HeatingEnergyTwoColourPt65, YieldTwoColourPt65./YieldTwoColourPt65(1))
scatter(HeatingEnergyArFPt5, YieldArFPt5./YieldArFPt5(1))
hold off
set(gca, 'YScale', 'log')

figure
plot(HeatingEnergyPt25./1000, YieldPt25./YieldPt25(1))
hold on
plot(HeatingEnergyPt35./1000, YieldPt35./YieldPt35(1))
plot(HeatingEnergyPt5./1000, YieldPt5./YieldPt5(1))
plot(HeatingEnergyPt65./1000, YieldPt65./YieldPt65(1))
plot(HeatingEnergyTwoColourPt5./1000, YieldTwoColourPt5./YieldTwoColourPt5(1))
plot(HeatingEnergyTwoColourPt65./1000, YieldTwoColourPt65./YieldTwoColourPt65(1))
plot(HeatingEnergyArFPt5./1000, YieldArFPt5./YieldArFPt5(1))
hold off
set(gca, 'YScale', 'log')
ylabel('Relative yield amplification')
xlabel('Heating energy (kJ)')

[X, Y1] = ArrayReturn(YieldPt25, YieldPt35, YieldPt5, YieldPt65, YieldTwoColourPt5, YieldTwoColourPt65, YieldArFPt5, 5);
[X, Y2] = ArrayReturn(YieldPt25, YieldPt35, YieldPt5, YieldPt65, YieldTwoColourPt5, YieldTwoColourPt65, YieldArFPt5, 15);
[X, Y3] = ArrayReturn(YieldPt25, YieldPt35, YieldPt5, YieldPt65, YieldTwoColourPt5, YieldTwoColourPt65, YieldArFPt5, 30);

figure
scatter(X, Y1)
hold on
scatter(X, Y2)
scatter(X, Y3)
hold off
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

%% Function
function [HeatingEnergySave, LaserEnergySave, GainSave, YieldSave, TotalInputEnergySave, GainWithHeatingSave, BurnUpSave] = MassLoad(file, ConversionAdjustment)
load(file);
HeatingEnergySave = HeatingEnergyAll;
LaserEnergySave = RequiredLaserEnergyAll;
GainSave = GainAll;
YieldSave = GainSave.*LaserEnergySave;
TotalInputEnergySave = LaserEnergySave + (ConversionAdjustment.*HeatingEnergySave.');
GainWithHeatingSave = YieldSave./TotalInputEnergySave;
BurnUpSave = BurnUpAll;
end

function [X, Y] = ArrayReturn(YieldPt25, YieldPt35, YieldPt5, YieldPt65, YieldTwoColourPt5, YieldTwoColourPt65, YieldArFPt5, Index)
X = [YieldPt25(1), YieldPt35(1), YieldPt5(1), YieldPt65(1), YieldTwoColourPt5(1), YieldTwoColourPt65(1), YieldArFPt5(1)];
Y = [YieldPt25(Index), YieldPt35(Index), YieldPt5(Index), YieldPt65(Index), YieldTwoColourPt5(Index), YieldTwoColourPt65(Index), YieldArFPt5(Index)];
Y = Y./X;
end