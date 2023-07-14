%Script to produce Hyades Input decks. Parameters are specified here, and
%additional scripts are called to do the meshing and write the input deck.

%Multiple configurations can be stored by putting into the if/else loop as
%'modes', and the 'mode' being specified in the first line.

%Options for sequences are:

%Change thicknesses - set the thicknesses as arrays. Each possible
%combination of the different arrays will be ran.

%Change laser timings (RelativeTimingMode OFF) - set the different laser
%timings in the arrays. All possible combinations will be ran. Arrays can
%be set as constants relative to each other.

%Change laser timings (RelativeTimingMode ON) - set the different laser
%timings in the arrays. No combinations will be performed: the arrays must
%be the same length, and a run will be performed using all the 1st values,
%then all the second values, etc. All the permutations of two arrays can be
%calculated using AbsolutePulseTimings. One array can be varied relative to
%another using RelativePulseTimings (and the first will be changed so that
%the arrays have the same lengths and pairings). If one of the arrays is a
%single constant value, this will be automatically replciated to the other
%array lengths.

%Other parameters (laser end times, vapour densities, electron heating,
%second laser timigns) can also be varied.

%FoamOnly makes the second layer a dry CD foam (no wetting). TNBurn enables
%the thermonuclear burn in layers 1 and 2. Electron Heating and Second
%Laser turn these two things on.

%SizeMod is used to ensure the laser power scales correctly with size.
%Rather than changing radius directly, scale it to a 0.285cm capsule
%(sizemod=1) and multiply by sizemod. Sizemod is then used to adjust the
%laser power accordingly. 
 
%Not all the old examples will work - may need to check that the final
%function contains all the variables, and that all the second laser
%variables are defined (even if they are just NaN).

%Mode is set up so that you can have multiple configurations saved in the
%script, and select which one you want to run
Mode=2.6;


%Mode 1 - Explanation of what info you need. Produces a single input deck
%called File1 in the MultiFile folder. If there is already a File1, it is
%overwritten.
if Mode==1
    
        %Wavelength - 0.351 for two colour, 0.193 for ArF
    Wavelength = 0.351;
    
    %SizeMod is required. 
     SizeMod = 0.75;
   %Relative TimingsMode=0 normally, or =1 if you have used RelativePulseTimings.
    RelativeTimingsMode = 1;
    %FoamOnly=1 means no DT in the foam layer. Leave as =0.
    FoamOnly=0;
    %TNBurn=1 activates the TN burn package (fusion reactions)
    TNBurn=1;
       
    %Set your pulse timings (i.e. switch on time of each pulse)
    Pulse1Time = 0;
    Pulse2Time = [3.6].*10^-9;
    Pulse3Time = [6.8].*10^-9;
    Pulse4Time = [8.2].*10^-9;

    %Pulse Powers are automatically set using sizemod so they satisfy our
    %intensity constraint.
    Pulse1Power = 2**((0.351/Wavelength)^2)*SizeMod^2;
    Pulse2Power = 14**((0.351/Wavelength)^2)*SizeMod^2;;
    Pulse3Power = 98**((0.351/Wavelength)^2)*SizeMod^2;
    Pulse4Power = 692.13**((0.351/Wavelength)^2)*SizeMod^2;
    
    %Other laser settings. We keep rise time (how long to increase power as
    %0.2ns. Laser End Time is when we turn the lasers off.
    RiseTime = 2e-10;
    LaserEndTime =[14].*10^-9 * (SizeMod/0.65);
    
    %Set the capsule dimensions (cm)
%     RadiusVapour = [0.165:0.0005:0.169] .* (SizeMod/0.65);
    RadiusVapour = [0.16705] .* (SizeMod/0.65);
%     RadiusIce = [0.181:0.00025:0.184] .* (SizeMod/0.65);
    RadiusIce = [0.182] .* (SizeMod/0.65);
    RadiusCH=0.18525 * (SizeMod/0.65);
    DensityVapour=[0.0009];
    
    %Set time spacing for the simulation. At ChangeTime(i), the TimeSpacing
    %is set to TimeSpacing(i). The final time in ChangeTimeAt is the end
    %time for the simulation.
    ChangeTimeAt = [0 15 24].*10^-9;
    TimeSpacing = [1E-10 1E-11];
    
    %ElectronHeatingOn=0 for no electron heating - all other values =NaN
    %unless ElectronHeatingOn=1
    ElectronHeatingOn=0;
    HeatingTimings = NaN;
    Heating = NaN;
    HeatingZones = NaN;
    HeatingEnergy = NaN;
        
    %SecondLaser=1 only if we are using the ArF laser. Like for electron
    %heating, other values are NaN unless SecondLaser=1
    SecondLaser=0;
    SecondLaserWavelength = NaN;
    SecondLaserOnTime = NaN;
    SecondLaserOffTime = NaN;
    [SecondLaserPower] = NaN;
    
    %Runs the function which will mesh the problem and write the input
    %deck. Does not require any changing, but must be included.
      [NumLayers, MassDiff, Method, Regions, ...
    Time, PulseOnTime, PulsePower] = GeneratorFunction(...
    RadiusVapour, RadiusIce, RadiusCH, DensityVapour, FoamOnly, TNBurn,...
    Pulse1Time, Pulse2Time, Pulse3Time, Pulse4Time, Pulse1Power, Pulse2Power, Pulse3Power, Pulse4Power, ...
    RelativeTimingsMode, RiseTime, LaserEndTime, ChangeTimeAt, TimeSpacing,...
    ElectronHeatingOn, HeatingTimings, Heating, HeatingZones, HeatingEnergy,...
     SecondLaser, SecondLaserOnTime, SecondLaserOffTime, SecondLaserWavelength, SecondLaserPower, Wavelength);
end


%Mode 2 - Demonstration of second laser, and of changing thickness.
if Mode==2
    
        %Wavelength - 0.351 for two colour, 0.193 for ArF
    Wavelength = 0.351;
    
    SizeMod = 0.5;
    RelativeTimingsMode = 0;
    FoamOnly=0;
    TNBurn=1;
      
    Pulse1Time = 0;
    Pulse2Time = [0.9].*10^-9;
    Pulse3Time = [3.8].*10^-9;
    Pulse4Time = [6.2].*10^-9;
    SecondLaserOnTime = [9.2]*10^-9;

%       Pulse1Time = 0;
%     Pulse3Time = [3.6:0.1:4.6].*10^-9;
%     [Pulse3Time, SecondLaserOnTime] = RelativePulseTiming(Pulse3Time, 4, 5, 0.1) ;  
%     Pulse2Time = Pulse3Time - 2.4*10^-9;
%     Pulse4Time = [5]*10^-9;

    Pulse1Power = 2**((0.351/Wavelength)^2)*SizeMod^2;
    Pulse2Power = 14**((0.351/Wavelength)^2)*SizeMod^2;;
    Pulse3Power = 98**((0.351/Wavelength)^2)*SizeMod^2;
    Pulse4Power = 692.13**((0.351/Wavelength)^2)*SizeMod^2;
    RiseTime = 2e-10;
    LaserEndTime = [11].*10^-9;
    
    %There are two arrays for thickness. An input deck will be produced for
    %every combination of these two arrays. (You can only vary one thing at
    %a time. While this run includes a second laser, it isn't changing
    %between runs - so the thickness can be changed and doesn't need to be
    %constant).
%     RadiusVapour = [0.123:0.0005:0.127];
%     RadiusIce = [0.135:0.0005:0.139];
    RadiusVapour = [0.125];
    RadiusIce = [0.137];
    RadiusCH=0.1425;
    DensityVapour=[0.0012];
    ChangeTimeAt = [0 5 8 11 15].*10^-9;
    TimeSpacing = [1E-10 1E-11 1E-10 1E-11];
    
    ElectronHeatingOn=0;
    HeatingTimings = NaN;
    Heating = NaN;
    HeatingZones = NaN;
    HeatingEnergy = NaN;
    
    %SecondLaser is now on. Specify the wavelength as a single value.
    %OnTime, OffTime and Radius must all have the same number of elements -
    %if there are i elements you will produce i input decks. File 1 will
    %use the first value from each array, File2 the second value, etc. You
    %can use repmat to repeat a single value so that you have the right
    %number of elements, if that value isn't changing over the input decks.
    %Radius should stay as RadiusCH (i.e. the outer radius). SecondLaserFn
    %determines the correct power, given the wavelength and radius.
    SecondLaser=1;
    SecondLaserWavelength = 0.193; %ArF frequency
    SecondLaserOffTime = LaserEndTime;
    SecondLaserRadius = [RadiusCH];
    [SecondLaserPower] = SecondLaserFn(SecondLaserWavelength, SecondLaserRadius);
  
    [NumLayers, MassDiff, Method, Regions, ...
    Time, PulseOnTime, PulsePower] = GeneratorFunction(...
    RadiusVapour, RadiusIce, RadiusCH, DensityVapour, FoamOnly, TNBurn,...
    Pulse1Time, Pulse2Time, Pulse3Time, Pulse4Time, Pulse1Power, Pulse2Power, Pulse3Power, Pulse4Power, ...
    RelativeTimingsMode, RiseTime, LaserEndTime, ChangeTimeAt, TimeSpacing,...
    ElectronHeatingOn, HeatingTimings, Heating, HeatingZones, HeatingEnergy,...
    SecondLaser, SecondLaserOnTime, LaserEndTime, SecondLaserWavelength, SecondLaserPower, Wavelength);
end



%Mode 2.5 - Demonstration of second laser, and of changing thickness.
if Mode==2.5
    
        %Wavelength - 0.351 for two colour, 0.193 for ArF
    Wavelength = 0.351;
    
    SizeMod = 0.5;
    RelativeTimingsMode = 0;
    FoamOnly=0;
    TNBurn=1;
    UnheatedFile = '\\aldaq1.physics.ox.ac.uk\Archer\Robert\Tat summer project\SummerProjectsMatlab\ElectronHeating_thirdharmonic\0.5Size4PulseCutoff.cdf';
    
    Pulse1Time = 0;
    Pulse2Time = [0.9].*10^-9;
    Pulse3Time = [3.8].*10^-9;
    Pulse4Time = [6.2].*10^-9;
    SecondLaserOnTime = [9.2]*10^-9;

%       Pulse1Time = 0;
%     Pulse3Time = [3.6:0.1:4.6].*10^-9;
%     [Pulse3Time, SecondLaserOnTime] = RelativePulseTiming(Pulse3Time, 4, 5, 0.1) ;  
%     Pulse2Time = Pulse3Time - 2.4*10^-9;
%     Pulse4Time = [5]*10^-9;

    Pulse1Power = 2**((0.351/Wavelength)^2)*SizeMod^2;
    Pulse2Power = 14**((0.351/Wavelength)^2)*SizeMod^2;;
    Pulse3Power = 98**((0.351/Wavelength)^2)*SizeMod^2;
    Pulse4Power = 692.13**((0.351/Wavelength)^2)*SizeMod^2;
    RiseTime = 2e-10;
    LaserEndTime = [11].*10^-9;
    
    %There are two arrays for thickness. An input deck will be produced for
    %every combination of these two arrays. (You can only vary one thing at
    %a time. While this run includes a second laser, it isn't changing
    %between runs - so the thickness can be changed and doesn't need to be
    %constant).
%     RadiusVapour = [0.123:0.0005:0.127];
%     RadiusIce = [0.135:0.0005:0.139];
    RadiusVapour = [0.125];
    RadiusIce = [0.137];
    RadiusCH=0.1425;
    DensityVapour=[0.0012];
    ChangeTimeAt = [0 5 8 11 15].*10^-9;
    TimeSpacing = [1E-10 1E-11 1E-10 1E-11];
    
    %Turns the electron heating on ( = 1) and off (=0)
    ElectronHeatingOn=1;
    HeatingTimings = [];
    Heating = [];
    HeatingZones = [];
    ChangeTimeAt = [];
    TimeSpacing = [];
    HeatingkJ = 4;
    HeatingEnergy = [];
    [BangTime, Time, TotalVapourMass, IceBoundaryIndex] = UnheatedFileReader(UnheatedFile)
    %HSMass = 8.7901E-5;
   
    %Calculates the appropriate amount of electron heating, and specifies
    %the different times that the heating should be applied.
    for HeatingStartTime = [(BangTime*10^9-0.4) :0.025: (BangTime*10^9+0.2)];
    [~, TimeIndex] = min(abs(Time - HeatingStartTime*10^-9)); %Total Vapour mass is actually the same regardless of time, but done this for continuity with above.
    HeatingConv = HeatingkJ*(10^10)/(5E-13 * TotalVapourMass(TimeIndex)); %EnergyKJ * KJ2Erg / (Time*Mass) to give Erg/g/sec
    HeatingTimings = [HeatingTimings; HeatingStartTime HeatingStartTime+0.0001 HeatingStartTime+0.0005 HeatingStartTime+0.0006];
    Heating = [Heating; 0 HeatingConv HeatingConv 0];
    HeatingZones = [HeatingZones; IceBoundaryIndex];
    ChangeTimeAt = [ChangeTimeAt; 0 11 HeatingStartTime HeatingStartTime+0.002 HeatingStartTime+0.02 Time(end).*10^9];
    TimeSpacing = [TimeSpacing;  1E-10 1E-11 1E-13 1E-12 2E-11 ];
    HeatingEnergy = [HeatingEnergy; HeatingkJ];
    end
    HeatingTimings=HeatingTimings*10^-9;
    ChangeTimeAt=ChangeTimeAt*10^-9;
    HeatingEnergy = HeatingEnergy*1000;
    
    %SecondLaser is now on. Specify the wavelength as a single value.
    %OnTime, OffTime and Radius must all have the same number of elements -
    %if there are i elements you will produce i input decks. File 1 will
    %use the first value from each array, File2 the second value, etc. You
    %can use repmat to repeat a single value so that you have the right
    %number of elements, if that value isn't changing over the input decks.
    %Radius should stay as RadiusCH (i.e. the outer radius). SecondLaserFn
    %determines the correct power, given the wavelength and radius.
    SecondLaser=1;
    SecondLaserWavelength = 0.193; %ArF frequency
    SecondLaserOffTime = LaserEndTime;
    SecondLaserRadius = [RadiusCH];
    [SecondLaserPower] = SecondLaserFn(SecondLaserWavelength, SecondLaserRadius);
  
    [NumLayers, MassDiff, Method, Regions, ...
    Time, PulseOnTime, PulsePower] = GeneratorFunction(...
    RadiusVapour, RadiusIce, RadiusCH, DensityVapour, FoamOnly, TNBurn,...
    Pulse1Time, Pulse2Time, Pulse3Time, Pulse4Time, Pulse1Power, Pulse2Power, Pulse3Power, Pulse4Power, ...
    RelativeTimingsMode, RiseTime, LaserEndTime, ChangeTimeAt, TimeSpacing,...
    ElectronHeatingOn, HeatingTimings, Heating, HeatingZones, HeatingEnergy,...
    SecondLaser, SecondLaserOnTime, LaserEndTime, SecondLaserWavelength, SecondLaserPower, Wavelength);
end

%Mode 2.5 - Demonstration of second laser, and of changing thickness.
if Mode==2.6
    
        %Wavelength - 0.351 for two colour, 0.193 for ArF
    Wavelength = 0.351;
    
    SizeMod = 0.5;
    RelativeTimingsMode = 0;
    FoamOnly=0;
    TNBurn=1;
    UnheatedFile = '\\aldaq1.physics.ox.ac.uk\Archer\Robert\Tat summer project\SummerProjectsMatlab\ElectronHeating_thirdharmonic\0.5Size4PulseCutoff.cdf';
    
    Pulse1Time = 0;
    Pulse2Time = [0.9].*10^-9;
    Pulse3Time = [3.8].*10^-9;
    Pulse4Time = [6.2].*10^-9;
    SecondLaserOnTime = [9.2]*10^-9;

%       Pulse1Time = 0;
%     Pulse3Time = [3.6:0.1:4.6].*10^-9;
%     [Pulse3Time, SecondLaserOnTime] = RelativePulseTiming(Pulse3Time, 4, 5, 0.1) ;  
%     Pulse2Time = Pulse3Time - 2.4*10^-9;
%     Pulse4Time = [5]*10^-9;

    Pulse1Power = 2**((0.351/Wavelength)^2)*SizeMod^2;
    Pulse2Power = 14**((0.351/Wavelength)^2)*SizeMod^2;;
    Pulse3Power = 98**((0.351/Wavelength)^2)*SizeMod^2;
    Pulse4Power = 692.13**((0.351/Wavelength)^2)*SizeMod^2;
    RiseTime = 2e-10;
    LaserEndTime = [11].*10^-9;
    
    %There are two arrays for thickness. An input deck will be produced for
    %every combination of these two arrays. (You can only vary one thing at
    %a time. While this run includes a second laser, it isn't changing
    %between runs - so the thickness can be changed and doesn't need to be
    %constant).
%     RadiusVapour = [0.123:0.0005:0.127];
%     RadiusIce = [0.135:0.0005:0.139];
    RadiusVapour = [0.125];
    RadiusIce = [0.137];
    RadiusCH=0.1425;
    DensityVapour=[0.0012];
    ChangeTimeAt = [0 5 8 11 15].*10^-9;
    TimeSpacing = [1E-10 1E-11 1E-10 1E-11];
    
    %Turns the electron heating on ( = 1) and off (=0)
    ElectronHeatingOn=1;
    HeatingTimings = [];
    Heating = [];
    HeatingZones = [];
    ChangeTimeAt = [];
    TimeSpacing = [];
    HeatingkJ = 4;
    HeatingEnergy = [];
    %HSMass = 8.7901E-5;
   
    %From the changing timing similations, you should find the optimal time
    %for the heating to occur. Put that in here.
    HeatingStartTime = 13.2552;
    
    %This code then generates the right variables so that the later code
    %will produce a series of files of varying power.
    for HeatingkJ = [2:2:60] %Original code
    %for HeatingkJ = [0:0.2:5]
    [~, TimeIndex] = min(abs(Time - HeatingStartTime*10^-9)); %Total Vapour mass is actually the same regardless of time, but done this for continuity with above.
    HeatingConv = HeatingkJ*(10^10)/(5E-13 * TotalVapourMass(TimeIndex)); %EnergyKJ * KJ2Erg / (Time*Mass) to give Erg/g/sec
    HeatingTimings = [HeatingTimings; HeatingStartTime HeatingStartTime+0.0001 HeatingStartTime+0.0005 HeatingStartTime+0.0006];
    Heating = [Heating; 0 HeatingConv HeatingConv 0];
    HeatingZones = [HeatingZones; IceBoundaryIndex];
    ChangeTimeAt = [ChangeTimeAt; 0 11 HeatingStartTime HeatingStartTime+0.002 HeatingStartTime+0.02 Time(end).*10^9];
    TimeSpacing = [TimeSpacing;  1E-10 1E-11 1E-13 1E-12 2E-11 ];
    HeatingEnergy = [HeatingEnergy; HeatingkJ];
    end
    HeatingTimings=HeatingTimings*10^-9;
    ChangeTimeAt=ChangeTimeAt*10^-9;
    HeatingEnergy = HeatingEnergy*1000;
    
    %SecondLaser is now on. Specify the wavelength as a single value.
    %OnTime, OffTime and Radius must all have the same number of elements -
    %if there are i elements you will produce i input decks. File 1 will
    %use the first value from each array, File2 the second value, etc. You
    %can use repmat to repeat a single value so that you have the right
    %number of elements, if that value isn't changing over the input decks.
    %Radius should stay as RadiusCH (i.e. the outer radius). SecondLaserFn
    %determines the correct power, given the wavelength and radius.
    SecondLaser=1;
    SecondLaserWavelength = 0.193; %ArF frequency
    SecondLaserOffTime = LaserEndTime;
    SecondLaserRadius = [RadiusCH];
    [SecondLaserPower] = SecondLaserFn(SecondLaserWavelength, SecondLaserRadius);
  
    [NumLayers, MassDiff, Method, Regions, ...
    Time, PulseOnTime, PulsePower] = GeneratorFunction(...
    RadiusVapour, RadiusIce, RadiusCH, DensityVapour, FoamOnly, TNBurn,...
    Pulse1Time, Pulse2Time, Pulse3Time, Pulse4Time, Pulse1Power, Pulse2Power, Pulse3Power, Pulse4Power, ...
    RelativeTimingsMode, RiseTime, LaserEndTime, ChangeTimeAt, TimeSpacing,...
    ElectronHeatingOn, HeatingTimings, Heating, HeatingZones, HeatingEnergy,...
    SecondLaser, SecondLaserOnTime, LaserEndTime, SecondLaserWavelength, SecondLaserPower, Wavelength);
end






%Mode 3 - Using relative pulse timings to change one pulse relative to another.
if Mode == 3
    
    SizeMod = 0.5;
    RelativeTimingsMode = 1;
    FoamOnly=0;
    TNBurn=1;
   
    
    %Laser timings are changed. Pulse 2 time is constant. An array is given
    %for pulse 3 times. Pulse 4 times are set relative to pulse 3 times -
    %for each pulse 3 value, input decks will be produced for pulse 4
    %values ranging from 0.7ns later to 1.5ns later, in 0.1ns intervals.
    %RelativeTimingsMode=1 must be on for this to work. I only ever vary
    %two sets of the timings in a given simulation run. We could have used
    %a single value from Pulse2Time (i.e. Pulse2Time = 2ns), but instead we
    %have set this relative to Pulse 3 Time.
       Pulse1Time = 0;
    Pulse3Time = [3.6:0.1:4.6].*10^-9;
    [Pulse3Time, Pulse4Time] = RelativePulseTiming(Pulse3Time, 0.7, 1.5, 0.1) ;  
    Pulse2Time = Pulse3Time - 2.4*10^-9;
 
    Pulse1Power = 2**((0.351/Wavelength)^2)*SizeMod^2;
    Pulse2Power = 14**((0.351/Wavelength)^2)*SizeMod^2;;
    Pulse3Power = 98**((0.351/Wavelength)^2)*SizeMod^2;
    Pulse4Power = 692.13**((0.351/Wavelength)^2)*SizeMod^2;
    RiseTime = 2e-10;
    LaserEndTime = 11.*10^-9;
    
    RadiusVapour = 0.125;
    RadiusIce = 0.14;
    RadiusCH=0.1425;
    DensityVapour=[0.0011];
    ChangeTimeAt = [0 6.5 8 11 16 20].*10^-9;
    TimeSpacing = [1E-10 1E-11 1E-10 1E-11 2E-11];
    
    ElectronHeatingOn=0;
    HeatingTimings = NaN;
    Heating = NaN;
    HeatingZones = NaN;
    HeatingEnergy = NaN;
    
    SecondLaser=0;
    SecondLaserWavelength = NaN;
    SecondLaserOnTime = NaN;
    SecondLaserOffTime = NaN;
    [SecondLaserPower] = NaN;
    
      
    [NumLayers, MassDiff, Method, Regions, ...
    Time, PulseOnTime, PulsePower] = GeneratorFunction(...
    RadiusVapour, RadiusIce, RadiusCH, DensityVapour, FoamOnly, TNBurn,...
    Pulse1Time, Pulse2Time, Pulse3Time, Pulse4Time, Pulse1Power, Pulse2Power, Pulse3Power, Pulse4Power, ...
    RelativeTimingsMode, RiseTime, LaserEndTime, ChangeTimeAt, TimeSpacing,...
    ElectronHeatingOn, HeatingTimings, Heating, HeatingZones, HeatingEnergy,...
    SecondLaser, SecondLaserOnTime, SecondLaserOffTime, SecondLaserWavelength, SecondLaserPower);
end


%Mode 4 - Vary pulse 2 and pulse 3 absolutely, no pulse 4
if Mode==4
    RelativeTimingsMode=0;
    FoamOnly=0;
    TNBurn=1;
    SizeMod = 0.35;
    
    %All combinations of pulse 2 time and pulse 3 time will be ran. By
    %setting pulse 4 time after the end of the simulation, pulse 4 will not
    %come on (and thus we only have 3 pulses)
    Pulse1Time = 0;
    Pulse2Time = [2.1:0.1:2.9].*10^-9;
    Pulse3Time = [3.4:0.1:5.4].*10^-9;
    Pulse4Time = [20].*10^-9;
    
    Pulse1Power = 2*SizeMod^2;
    Pulse2Power = 14*SizeMod^2;
    Pulse3Power = 98*SizeMod^2;
    Pulse4Power = 692.13*SizeMod^2;
    RiseTime = 2e-10;
    LaserEndTime = 35.*10^-9;
    RadiusVapour = [0.251].*SizeMod;
    RadiusIce = [0.278].*SizeMod;
    RadiusCH=0.285*SizeMod;
    DensityVapour=0.001;
    ChangeTimeAt = [0 3 5 10]*10^-9;
    TimeSpacing = [1E-11 2E-11 1E-11];
    
    ElectronHeatingOn=0;
    HeatingTimings = NaN;
    Heating = NaN;
    HeatingZones = NaN;
    HeatingEnergy = NaN;
    
    SecondLaser=0;
    SecondLaserWavelength = NaN;
    SecondLaserOnTime = NaN;
    SecondLaserOffTime = NaN;
    [SecondLaserPower] = NaN;
    
      
    [NumLayers, MassDiff, Method, Regions, ...
    Time, PulseOnTime, PulsePower] = GeneratorFunction(...
    RadiusVapour, RadiusIce, RadiusCH, DensityVapour, FoamOnly, TNBurn,...
    Pulse1Time, Pulse2Time, Pulse3Time, Pulse4Time, Pulse1Power, Pulse2Power, Pulse3Power, Pulse4Power, ...
    RelativeTimingsMode, RiseTime, LaserEndTime, ChangeTimeAt, TimeSpacing,...
    ElectronHeatingOn, HeatingTimings, Heating, HeatingZones, HeatingEnergy,...
    SecondLaser, SecondLaserOnTime, SecondLaserOffTime, SecondLaserWavelength, SecondLaserPower);
end


%Mode 5 - Recreate Tabak paper
if Mode==5
    
    %SizeMod is required. 
     SizeMod = 0.1410/0.285;
    %Relative TimingsMode=0 normally, or =1 if you have used RelativePulseTimings.
    RelativeTimingsMode = 0;
    %FoamOnly=1 means no DT in the foam layer. Leave as =0.
    FoamOnly=0;
    %TNBurn=1 activates the TN burn package (fusion reactions)
    TNBurn=1;
   
    
    %Set your pulse timings (i.e. switch on time of each pulse)
    Pulse1Time = 0;
    Pulse2Time = [2.2].*10^-9;
    Pulse3Time = [4.6].*10^-9;
    %Pulse4Time = [5.5].*10^-9;
    Pulse4Time = NaN;

    %Pulse Powers are automatically set using sizemod so they satisfy our
    %intensity constraint.
%     Pulse1Power = 2*SizeMod^2;
%     Pulse2Power = 14*SizeMod^2;
%     Pulse3Power = 98*SizeMod^2;
%     Pulse4Power = 692.13*SizeMod^2;
Pulse1Power = 2E12*4*pi()*(0.1410^2)/10^12;
Pulse2Power = 2E13*4*pi()*(0.1410^2)/10^12;
Pulse3Power = 2.5E15*4*pi()*(0.1410^2)/10^12;
Pulse4Power=NaN;
    
    %Other laser settings. We keep rise time (how long to increase power as
    %0.2ns. Laser End Time is when we turn the lasers off.
    RiseTime = 2e-10;
    LaserEndTime =35.*10^-9;
    
    %Set the capsule dimensions (cm)
    RadiusVapour = 0.1;
    RadiusIce = 0.123;
    RadiusCH=0.1410;
    DensityVapour=[0.0135];
    
    %Set time spacing for the simulation. At ChangeTime(i), the TimeSpacing
    %is set to TimeSpacing(i). The final time in ChangeTimeAt is the end
    %time for the simulation.
    ChangeTimeAt = [0 50].*10^-9;
    TimeSpacing = [1E-10];
    
    %ElectronHeatingOn=0 for no electron heating - all other values =NaN
    %unless ElectronHeatingOn=1
    ElectronHeatingOn=0;
    HeatingTimings = NaN;
    Heating = NaN;
    HeatingZones = NaN;
    HeatingEnergy = NaN;
    
    
    %SecondLaser=1 only if we are using the ArF laser. Like for electron
    %heating, other values are NaN unless SecondLaser=1
    SecondLaser=0;
    SecondLaserWavelength = NaN;
    SecondLaserOnTime = NaN;
    SecondLaserOffTime = NaN;
    [SecondLaserPower] = NaN;
    
    %Runs the function which will mesh the problem and write the input
    %deck. Does not require any changing, but must be included.
      [NumLayers, MassDiff, Method, Regions, ...
    Time, PulseOnTime, PulsePower] = GeneratorFunction(...
    RadiusVapour, RadiusIce, RadiusCH, DensityVapour, FoamOnly, TNBurn,...
    Pulse1Time, Pulse2Time, Pulse3Time, Pulse4Time, Pulse1Power, Pulse2Power, Pulse3Power, Pulse4Power, ...
    RelativeTimingsMode, RiseTime, LaserEndTime, ChangeTimeAt, TimeSpacing,...
    ElectronHeatingOn, HeatingTimings, Heating, HeatingZones, HeatingEnergy,...
     SecondLaser, SecondLaserOnTime, SecondLaserOffTime, SecondLaserWavelength, SecondLaserPower);
end


%Mode 6 - Aligned Shock timings
if Mode==6
    RelativeTimingsMode=0;
    FoamOnly=0;
    TNBurn=1;
    SizeMod = 0.5;
    
    %All combinations of pulse 2 time and pulse 3 time will be ran. By
    %setting pulse 4 time after the end of the simulation, pulse 4 will not
    %come on (and thus we only have 3 pulses)
    Pulse1Time = 0;
    Pulse2Time = [7.615].*10^-9;
    Pulse3Time = [9.42].*10^-9;
    Pulse4Time = [10.05].*10^-9;
    
    Pulse1Power = 2*SizeMod^2;
    Pulse2Power = 14*SizeMod^2;
    Pulse3Power = 98*SizeMod^2;
    Pulse4Power = 692.13*SizeMod^2;
    RiseTime = 2e-10;
    LaserEndTime = 18.*10^-9;
%     RadiusVapour = [0.251].*SizeMod;
    RadiusVapour = 0.1305;
%     RadiusIce = [0.278].*SizeMod;
    RadiusIce = 0.1395;
%     RadiusCH=0.285*SizeMod;
    RadiusCH=0.1425;

    DensityVapour=[0.0001:0.0001:0.001];
    ChangeTimeAt = [0 10 10.4 10.55 17]*10^-9;
    TimeSpacing = [1E-10 1E-11 1E-12 1E-11 ];
    
    ElectronHeatingOn=0;
    HeatingTimings = NaN;
    Heating = NaN;
    HeatingZones = NaN;
    HeatingEnergy = NaN;
    
    SecondLaser=0;
    SecondLaserWavelength = NaN;
    SecondLaserOnTime = NaN;
    SecondLaserOffTime = NaN;
    [SecondLaserPower] = NaN;
    
      
    [NumLayers, MassDiff, Method, Regions, ...
    Time, PulseOnTime, PulsePower] = GeneratorFunction(...
    RadiusVapour, RadiusIce, RadiusCH, DensityVapour, FoamOnly, TNBurn,...
    Pulse1Time, Pulse2Time, Pulse3Time, Pulse4Time, Pulse1Power, Pulse2Power, Pulse3Power, Pulse4Power, ...
    RelativeTimingsMode, RiseTime, LaserEndTime, ChangeTimeAt, TimeSpacing,...
    ElectronHeatingOn, HeatingTimings, Heating, HeatingZones, HeatingEnergy,...
    SecondLaser, SecondLaserOnTime, SecondLaserOffTime, SecondLaserWavelength, SecondLaserPower);
end

%Mode 7 - Recreate ArF implosion
if Mode==7
    RelativeTimingsMode=0;
    FoamOnly=0;
    TNBurn=1;
    SizeMod = 0.65;
    
    %All combinations of pulse 2 time and pulse 3 time will be ran. By
    %setting pulse 4 time after the end of the simulation, pulse 4 will not
    %come on (and thus we only have 3 pulses)
    Pulse1Time = 0;
    Pulse2Time = [7.615].*10^-9;
    Pulse3Time = [9.42].*10^-9;
    Pulse4Time = [10.05].*10^-9;
    
    Pulse1Power = 1.6538;
    Pulse2Power = 11.576;
    Pulse3Power = 81.037;
    Pulse4Power = 572.25;
    RiseTime = 2e-10;
    LaserEndTime = 18.*10^-9;
%     RadiusVapour = [0.251].*SizeMod;
    RadiusVapour = 0.1645;
%     RadiusIce = [0.278].*SizeMod;
    RadiusIce = 0.1818;
%     RadiusCH=0.285*SizeMod;
    RadiusCH=0.18525;

    DensityVapour=[0.0001:0.0001:0.001];
    ChangeTimeAt = [0 10 10.4 10.55 17]*10^-9;
    TimeSpacing = [1E-10 1E-11 1E-12 1E-11 ];
    
    ElectronHeatingOn=0;
    HeatingTimings = NaN;
    Heating = NaN;
    HeatingZones = NaN;
    HeatingEnergy = NaN;
    
    SecondLaser=0;
    SecondLaserWavelength = NaN;
    SecondLaserOnTime = NaN;
    SecondLaserOffTime = NaN;
    [SecondLaserPower] = NaN;
    
      
    [NumLayers, MassDiff, Method, Regions, ...
    Time, PulseOnTime, PulsePower] = GeneratorFunction(...
    RadiusVapour, RadiusIce, RadiusCH, DensityVapour, FoamOnly, TNBurn,...
    Pulse1Time, Pulse2Time, Pulse3Time, Pulse4Time, Pulse1Power, Pulse2Power, Pulse3Power, Pulse4Power, ...
    RelativeTimingsMode, RiseTime, LaserEndTime, ChangeTimeAt, TimeSpacing,...
    ElectronHeatingOn, HeatingTimings, Heating, HeatingZones, HeatingEnergy,...
    SecondLaser, SecondLaserOnTime, SecondLaserOffTime, SecondLaserWavelength, SecondLaserPower);
end

%Mode 8 - Recreate Connor Williams' capsule.
if Mode==8
    
    %SizeMod is required. 
     SizeMod = 0.178;
   %Relative TimingsMode=0 normally, or =1 if you have used RelativePulseTimings.
    RelativeTimingsMode = 0;
    %FoamOnly=1 means no DT in the foam layer. Leave as =0.
    FoamOnly=0;
    %TNBurn=1 activates the TN burn package (fusion reactions)
    TNBurn=1;
   
    
    %Set your pulse timings (i.e. switch on time of each pulse)
    Pulse1Time = 0;
    Pulse2Time = [0.4].*10^-9;
    Pulse3Time = [0.7].*10^-9;
    Pulse4Time = [1.0].*10^-9;

    %Pulse Powers are automatically set using sizemod so they satisfy our
    %intensity constraint.
    Pulse1Power = 27;
    Pulse2Power = 27;
    Pulse3Power = 27;
    Pulse4Power = 27;
    
    %Other laser settings. We keep rise time (how long to increase power as
    %0.2ns. Laser End Time is when we turn the lasers off.
    RiseTime = 2e-10;
    LaserEndTime =1.3.*10^-9;
    
    %Set the capsule dimensions (cm)
    RadiusVapour = 0.04665;
    RadiusIce = 0.05005;
    RadiusCH=0.0508;
    DensityVapour=[0.0003];
    
    %Set time spacing for the simulation. At ChangeTime(i), the TimeSpacing
    %is set to TimeSpacing(i). The final time in ChangeTimeAt is the end
    %time for the simulation.
    ChangeTimeAt = [0 3].*10^-9;
    TimeSpacing = [1E-11];
    
    %ElectronHeatingOn=0 for no electron heating - all other values =NaN
    %unless ElectronHeatingOn=1
    ElectronHeatingOn=0;
    HeatingTimings = NaN;
    Heating = NaN;
    HeatingZones = NaN;
    HeatingEnergy = NaN;
    
    
    %SecondLaser=1 only if we are using the ArF laser. Like for electron
    %heating, other values are NaN unless SecondLaser=1
    SecondLaser=0;
    SecondLaserWavelength = NaN;
    SecondLaserOnTime = NaN;
    SecondLaserOffTime = NaN;
    [SecondLaserPower] = NaN;
    
    %Runs the function which will mesh the problem and write the input
    %deck. Does not require any changing, but must be included.
      [NumLayers, MassDiff, Method, Regions, ...
    Time, PulseOnTime, PulsePower] = GeneratorFunction(...
    RadiusVapour, RadiusIce, RadiusCH, DensityVapour, FoamOnly, TNBurn,...
    Pulse1Time, Pulse2Time, Pulse3Time, Pulse4Time, Pulse1Power, Pulse2Power, Pulse3Power, Pulse4Power, ...
    RelativeTimingsMode, RiseTime, LaserEndTime, ChangeTimeAt, TimeSpacing,...
    ElectronHeatingOn, HeatingTimings, Heating, HeatingZones, HeatingEnergy,...
     SecondLaser, SecondLaserOnTime, SecondLaserOffTime, SecondLaserWavelength, SecondLaserPower);
end


%Mode 1 - Explanation of what info you need. Produces a single input deck
%called File1 in the MultiFile folder. If there is already a File1, it is
%overwritten.
if Mode==9
    
    %SizeMod is required. 
     SizeMod = 0.35;
   %Relative TimingsMode=0 normally, or =1 if you have used RelativePulseTimings.
    RelativeTimingsMode = 0;
    %FoamOnly=1 means no DT in the foam layer. Leave as =0.
    FoamOnly=1;
    %TNBurn=1 activates the TN burn package (fusion reactions)
    TNBurn=1;
   
    
    %Set your pulse timings (i.e. switch on time of each pulse)
    Pulse1Time = 0;
%     Pulse2Time = [3.6].*10^-9 * (SizeMod/0.65); %Original
%     Pulse3Time = [7.8].*10^-9* (SizeMod/0.65); %Original
%     Pulse3Time = [7:0.2:9].*10^-9* (SizeMod/0.65);
%     Pulse4Time = [9.5].*10^-9* (SizeMod/0.65); %Original
%     Pulse4Time = [9:0.2:11].*10^-9* (SizeMod/0.65);
Pulse2Time = [2.7].*10^-9;
Pulse3Time = [4.9].*10^-9;
Pulse4Time = [5.8].*10^-9;
% Pulse2Time = [3.5:0.2:5.1].*10^-9;
% Pulse3Time = [8:0.2:10].*10^-9;
% Pulse4Time = [10.9].*10^-9;


    %Pulse Powers are automatically set using sizemod so they satisfy our
    %intensity constraint.
    Pulse1Power = 2*SizeMod^2;
    Pulse2Power = 14*SizeMod^2;
    Pulse3Power = 98*SizeMod^2;
    Pulse4Power = 692.13*SizeMod^2;
    
    %Other laser settings. We keep rise time (how long to increase power as
    %0.2ns. Laser End Time is when we turn the lasers off.
    RiseTime = 2e-10;
    LaserEndTime =[8.5].*10^-9 ;
    
    %Set the capsule dimensions (cm)
%     RadiusVapour = [0.165:0.0005:0.169] .* (SizeMod/0.65);
    RadiusVapour = [0.0905];
%     RadiusIce = [0.181:0.00025:0.184] .* (SizeMod/0.65);
    RadiusIce = [0.098];
    RadiusCH=0.09975;
    DensityVapour=[0.00395:0.00001:0.00403];
    
    %Set time spacing for the simulation. At ChangeTime(i), the TimeSpacing
    %is set to TimeSpacing(i). The final time in ChangeTimeAt is the end
    %time for the simulation.
    ChangeTimeAt = [0 3 5 10].*10^-9;
    TimeSpacing = [1E-10 2E-11 1E-11];
    
    %ElectronHeatingOn=0 for no electron heating - all other values =NaN
    %unless ElectronHeatingOn=1
    ElectronHeatingOn=0;
    HeatingTimings = NaN;
    Heating = NaN;
    HeatingZones = NaN;
    HeatingEnergy = NaN;
    
    
    %SecondLaser=1 only if we are using the ArF laser. Like for electron
    %heating, other values are NaN unless SecondLaser=1
    SecondLaser=0;
    SecondLaserWavelength = NaN;
    SecondLaserOnTime = NaN;
    SecondLaserOffTime = NaN;
    [SecondLaserPower] = NaN;
    
    %Runs the function which will mesh the problem and write the input
    %deck. Does not require any changing, but must be included.
      [NumLayers, MassDiff, Method, Regions, ...
    Time, PulseOnTime, PulsePower] = GeneratorFunction(...
    RadiusVapour, RadiusIce, RadiusCH, DensityVapour, FoamOnly, TNBurn,...
    Pulse1Time, Pulse2Time, Pulse3Time, Pulse4Time, Pulse1Power, Pulse2Power, Pulse3Power, Pulse4Power, ...
    RelativeTimingsMode, RiseTime, LaserEndTime, ChangeTimeAt, TimeSpacing,...
    ElectronHeatingOn, HeatingTimings, Heating, HeatingZones, HeatingEnergy,...
     SecondLaser, SecondLaserOnTime, SecondLaserOffTime, SecondLaserWavelength, SecondLaserPower);
end

%Mode 1 - Explanation of what info you need. Produces a single input deck
%called File1 in the MultiFile folder. If there is already a File1, it is
%overwritten.
if Mode==10
    
    %SizeMod is required. 
     SizeMod = 1;
   %Relative TimingsMode=0 normally, or =1 if you have used RelativePulseTimings.
    RelativeTimingsMode = 0;
    %FoamOnly=1 means no DT in the foam layer. Leave as =0.
    FoamOnly=0;
    %TNBurn=1 activates the TN burn package (fusion reactions)
    TNBurn=1;
   
    
    %Set your pulse timings (i.e. switch on time of each pulse)
    Pulse1Time = 0;
%     Pulse2Time = [3.6].*10^-9 * (SizeMod/0.65); %Original
%     Pulse3Time = [7.8].*10^-9* (SizeMod/0.65); %Original
%     Pulse3Time = [7:0.2:9].*10^-9* (SizeMod/0.65);
%     Pulse4Time = [9.5].*10^-9* (SizeMod/0.65); %Original
%     Pulse4Time = [9:0.2:11].*10^-9* (SizeMod/0.65);
Pulse2Time = [4.2].*10^-9;
Pulse3Time = [8.8].*10^-9;
Pulse4Time = [10.9].*10^-9;
% Pulse2Time = [3.5:0.2:5.1].*10^-9;
% Pulse3Time = [8:0.2:10].*10^-9;
% Pulse4Time = [10.9].*10^-9;


    %Pulse Powers are automatically set using sizemod so they satisfy our
    %intensity constraint.
    Pulse1Power = 23;
    Pulse2Power = 230;
    Pulse3Power = 230;
    Pulse4Power = 230;
    
    %Other laser settings. We keep rise time (how long to increase power as
    %0.2ns. Laser End Time is when we turn the lasers off.
    RiseTime = 2e-10;
    LaserEndTime =[13].*10^-9;
    
    %Set the capsule dimensions (cm)
%     RadiusVapour = [0.165:0.0005:0.169] .* (SizeMod/0.65);
    RadiusVapour = [0.2008];
%     RadiusIce = [0.181:0.00025:0.184] .* (SizeMod/0.65);
    RadiusIce = [0.2338];
    RadiusCH=0.2343;
    DensityVapour=[0.001:0.0002:0.003];
    
    %Set time spacing for the simulation. At ChangeTime(i), the TimeSpacing
    %is set to TimeSpacing(i). The final time in ChangeTimeAt is the end
    %time for the simulation.
    ChangeTimeAt = [0 10 18].*10^-9;
    TimeSpacing = [1E-10 1E-11];
    
    %ElectronHeatingOn=0 for no electron heating - all other values =NaN
    %unless ElectronHeatingOn=1
    ElectronHeatingOn=0;
    HeatingTimings = NaN;
    Heating = NaN;
    HeatingZones = NaN;
    HeatingEnergy = NaN;
    
    
    %SecondLaser=1 only if we are using the ArF laser. Like for electron
    %heating, other values are NaN unless SecondLaser=1
    SecondLaser=0;
    SecondLaserWavelength = NaN;
    SecondLaserOnTime = NaN;
    SecondLaserOffTime = NaN;
    [SecondLaserPower] = NaN;
    
    %Runs the function which will mesh the problem and write the input
    %deck. Does not require any changing, but must be included.
      [NumLayers, MassDiff, Method, Regions, ...
    Time, PulseOnTime, PulsePower] = GeneratorFunction(...
    RadiusVapour, RadiusIce, RadiusCH, DensityVapour, FoamOnly, TNBurn,...
    Pulse1Time, Pulse2Time, Pulse3Time, Pulse4Time, Pulse1Power, Pulse2Power, Pulse3Power, Pulse4Power, ...
    RelativeTimingsMode, RiseTime, LaserEndTime, ChangeTimeAt, TimeSpacing,...
    ElectronHeatingOn, HeatingTimings, Heating, HeatingZones, HeatingEnergy,...
     SecondLaser, SecondLaserOnTime, SecondLaserOffTime, SecondLaserWavelength, SecondLaserPower);
end



%Mode 27 - Changing timings of 0.25 size capsule
if Mode==11
      %SizeMod is required. 
     SizeMod = 0.65;
   %Relative TimingsMode=0 normally, or =1 if you have used RelativePulseTimings.
    RelativeTimingsMode = 0;
    %FoamOnly=1 means no DT in the foam layer. Leave as =0.
    FoamOnly=0;
    %TNBurn=1 activates the TN burn package (fusion reactions)
    TNBurn=1;
    
    SizeMod = 0.65;
    Pulse1Time = 0;
    Pulse2Time = [3.6].*10^-9;
    Pulse3Time = [7.8].*10^-9;
%     [Pulse2Time, Pulse3Time] = RelativePulseTiming(Pulse2Time, 2, 3, 0.1);
    Pulse4Time = [9.5].*10^-9;
    %[Pulse3Time, Pulse4Time] = RelativePulseTiming(Pulse3Time, 0.7, 1.5, 0.1) ;  
    %Pulse4Time = 5.*10^-9

    Pulse1Power = 2*SizeMod^2;
    Pulse2Power = 14*SizeMod^2;
    Pulse3Power = 98*SizeMod^2;
    Pulse4Power = 692.13*SizeMod^2;
    RiseTime = 2e-10;
    LaserEndTime = [15].*10^-9;
    RadiusVapour = [0.16705];
    RadiusIce = [0.182];
    RadiusCH=0.285*SizeMod;
    DensityVapour=[0.001:0.0001:0.002];
  ChangeTimeAt = [0 8.5 11 13.5 17 20]*10^-9;
    TimeSpacing = [1E-10 1E-11 1E-10 1E-11 1E-10];

       
    %ElectronHeatingOn=0 for no electron heating - all other values =NaN
    %unless ElectronHeatingOn=1
    ElectronHeatingOn=0;
    HeatingTimings = NaN;
    Heating = NaN;
    HeatingZones = NaN;
    HeatingEnergy = NaN;
    
    
    %SecondLaser=1 only if we are using the ArF laser. Like for electron
    %heating, other values are NaN unless SecondLaser=1
    SecondLaser=0;
    SecondLaserWavelength = NaN;
    SecondLaserOnTime = NaN;
    SecondLaserOffTime = NaN;
    [SecondLaserPower] = NaN;
    
    %Runs the function which will mesh the problem and write the input
    %deck. Does not require any changing, but must be included.
      [NumLayers, MassDiff, Method, Regions, ...
    Time, PulseOnTime, PulsePower] = GeneratorFunction(...
    RadiusVapour, RadiusIce, RadiusCH, DensityVapour, FoamOnly, TNBurn,...
    Pulse1Time, Pulse2Time, Pulse3Time, Pulse4Time, Pulse1Power, Pulse2Power, Pulse3Power, Pulse4Power, ...
    RelativeTimingsMode, RiseTime, LaserEndTime, ChangeTimeAt, TimeSpacing,...
    ElectronHeatingOn, HeatingTimings, Heating, HeatingZones, HeatingEnergy,...
     SecondLaser, SecondLaserOnTime, SecondLaserOffTime, SecondLaserWavelength, SecondLaserPower);
end



function [BangTime, Time, TotalVapourMass, IceBoundaryIndex] = UnheatedFileReader(file)
%This code opens that previous file, so that some of the key settings can
%be imported. This section can be skipped!
Radius  = ncread(file,'R');
Time  = ncread(file,'DumpTimes');
Radiuscm  = ncread(file,'Rcm');
Density  = ncread(file,'Rho'); 
ElecTemp  = ncread(file,'Te')*1000; 
IonTemp  = ncread(file,'Ti')*1000; 
Neutrons = sum((0.8/(14.1 * 1.60218e-6))*ncread(file, 'Bpeprd'));
Neutrons(isnan(Neutrons)) = 0;
NeutronRate = sum((0.8/(14.1 * 1.60218e-6))*ncread(file, 'Bpeprdr'));
OldRhoR = Density.*Radius(1:end-1, :); 
Volume =  ncread(file,'Vol'); 
Mass = Volume(2:end,:).*Density;
Pressure = ncread(file,'Pres').*0.0000001; 
PressureGPa = Pressure.*10^-3;
TNOutput = ncread(file, 'Bpeprdr'); 
TNInput = ncread(file, 'Bpedep').*10^-7; 
TNInputRate = [zeros(size(TNInput,1),1),diff(TNInput,1,2)]; 
Velocity = ncread(file, 'U')./100000; 
IonThermalEnergy  = ncread(file,'Eion').*10^-7; 
ElectronThermalEnergy  = ncread(file,'Eelc').*10^-7; 
DepositedLaserPowerZones = ncread(file,'Deplas').*10^-7; 
DepositedLaserPower = sum(DepositedLaserPowerZones).';
DepositedEnergy = trapz(Time, DepositedLaserPower); 
LaserEnergy  = ncread(file,'Elasin')*(10^-7)/0.8;
LaserPower=[0; diff(LaserEnergy)./diff(Time)];
SimulationPower = LaserPower*0.8;
TotalLaserEnergy = LaserEnergy(end); 
TotalLaserPower = LaserPower;
LaserPowerDiff=(([0; diff(LaserPower)]./LaserPower)<0.0001).';
PulseIndices = findstr([0 LaserPowerDiff], [0 1])-1;
PulsePowers = LaserPower(PulseIndices);
m = (LaserPower(PulseIndices-1) - LaserPower(PulseIndices-2))./(Time(PulseIndices-1) - Time(PulseIndices-2));
c = LaserPower(PulseIndices-1) - m.*Time(PulseIndices-1);
PulseOnTime = ((LaserPower([1 PulseIndices(1:end-1)])-c)./m) - 0.0005E-7;
PulseOnAndRiseTime = (LaserPower([PulseIndices(1:end)])-c(1:end))./m(1:end)- 0.0005E-7;
Type=1;
Plots=0;
AnalysisCondensed
end


function [SecondLaserPower] = SecondLaserFn(SecondLaserWavelength, SecondLaserRadius)
 ILambda2 = 4E13;      
%ILambda2 = 8.307E13;
%ILambda2 = 1.6E14;
    SecondLaserIntensity = ILambda2./(SecondLaserWavelength.^2);
    SecondLaserPower = SecondLaserIntensity .* 4*pi() .* SecondLaserRadius.^2 / 10^12;
end


%(The absolutetiming function creates every possible permutation - if array 1 is
%[A B] and array 2 is [C D], this will give you all the possible
%combinations between 1 and 2 (i.e. A&C, A&D, B&C, B&D).)
function [PulseATime, PulseBTime] = AbsolutePulseTiming(PulseATimeValues, PulseBTimeValues)
    PulseATime=[];
    PulseBTime=[];
    for PulseACurrentValue=PulseATimeValues
        RollingPulseBTimeValues = PulseBTimeValues;
        RollingPulseATime = repmat(PulseACurrentValue, [1 length(RollingPulseBTimeValues)]);
        PulseATime = [PulseATime, RollingPulseATime];
        PulseBTime = [PulseBTime, RollingPulseBTimeValues];
    end
end


%(The relativetiming function allows you to create one array relative to
%another. If array 1 has values [A B], you can specify a setup where array
%2 can take a range of values relative to A (i.e. take the value in A and
%add everything from 1 to 10 in 1 increments for each value), and adjust A
%so they are paired up correctly).
function [PulseATime, PulseBTime] = RelativePulseTiming(PulseATimeValues, Lower, Upper, Increment)
    PulseATime=[];
    PulseBTime=[];
    Lower = Lower*10^-9;
    Upper = Upper*10^-9;
    Increment = Increment*10^-9;
    for PulseACurrentValue=PulseATimeValues
        RollingPulseBTimeValues = [PulseACurrentValue+Lower:Increment:PulseACurrentValue+Upper];
        RollingPulseATime = repmat(PulseACurrentValue, [1 length(RollingPulseBTimeValues)]);
        PulseATime = [PulseATime, RollingPulseATime];
        PulseBTime = [PulseBTime, RollingPulseBTimeValues];
    end
end

function [RadiusVapour, RadiusIce, RadiusCH] = AbsoluteThicknesses(RadiusVapourValues, RadiusIceValues, RadiusCHValues)
% Replicate the Radii to cover the possibiliies
%     RadiusVapour = repmat(RadiusVapourValues, [1 length(RadiusIceValues)])';
%     RadiusCH = repmat(RadiusCHValues, length(RadiusVapour), 1);
%     RadiusIce=repmat(RadiusIceValues,1,length(RadiusVapourValues))';
%     RadiusIce=RadiusIce(:)';

    %Replicate the Radii to cover the possibiliies. This is the simplest
    %way of doing it.
    RadiusVapour = repmat(RadiusVapourValues, [1 length(RadiusIceValues)]);
    RadiusCH = repmat(RadiusCHValues, length(RadiusVapour), 1);
    RadiusIce=repmat(RadiusIceValues',1,length(RadiusVapourValues))';
    RadiusIce=RadiusIce(:)';

end


function [NumLayers, MaxMassDiff, Method, Regions, ...
    Time, PulseOnTime, PulsePower] = GeneratorFunction(...
    RadiusVapour, RadiusIce, RadiusCH, DensityVapour, FoamOnly, TNBurn,...
    Pulse1Time, Pulse2Time, Pulse3Time, Pulse4Time, Pulse1Power, Pulse2Power, Pulse3Power, Pulse4Power, ...
    RelativeTimingsMode, RiseTime, LaserEndTime, ChangeTimeAt, TimeSpacing,...
    ElectronHeatingOn, HeatingTimings, Heating, HeatingZones, HeatingEnergy, ...
     SecondLaser, SecondLaserOnTime, SecondLaserOffTime, SecondLaserWavelength, SecondLaserPower, Wavelength)


%Mode1 - Thickness Scan
    if (length(RadiusIce)>1 || length(RadiusVapour)>1)
        disp("Varying Radius Ice and Radius Vapour")
        %Generate Radius Values
        [RadiusVapour, RadiusIce, RadiusCH] = AbsoluteThicknesses(RadiusVapour, RadiusIce, RadiusCH);
        %For each thickness
        for i=1:length(RadiusIce)
            %Mesh
             [NumLayers(i), MaxMassDiff(i), Method, Regions, ...
                 RadiusCHEndZone, LayersVapour, LayersIce, LayersCH, LayersCH1, LayersCHEndZone,...
                 RatioVapour, RatioIce, RatioCH, RatioCH1, RatioCHEndZone, ~, ~]...
                 = MeshSolver(RadiusVapour(i), RadiusIce(i), RadiusCH(i), DensityVapour, FoamOnly);
             %Write File
             FileWriter(RadiusVapour(i), RadiusIce(i), RadiusCH(i), RadiusCHEndZone, LayersVapour, LayersIce, LayersCH, LayersCH1, LayersCHEndZone, RatioVapour, RatioIce, RatioCH, RatioCH1, RatioCHEndZone, DensityVapour, i, Regions, FoamOnly,TNBurn,...
                  Pulse1Time, Pulse2Time, Pulse3Time, Pulse4Time, Pulse1Power,Pulse2Power,Pulse3Power, Pulse4Power, RiseTime, LaserEndTime, ChangeTimeAt, TimeSpacing, ElectronHeatingOn, HeatingTimings, Heating, HeatingZones, HeatingEnergy, ...
     SecondLaser, SecondLaserOnTime, SecondLaserOffTime, SecondLaserWavelength, SecondLaserPower, Wavelength);
              disp(['Meshing for file ', num2str(i), ' complete, number of layers = ', num2str(NumLayers(i)), ', max mass diff = ', num2str(MaxMassDiff(i))]);
        end
        figure
         plot(RadiusIce)
         hold on
        plot(RadiusVapour)
        plot(RadiusCH)
        hold off
        figure
        yyaxis left
        plot(NumLayers)
        yyaxis right
        plot(MaxMassDiff)
        if isnan(Pulse4Time)==1
        PulseOnTime = repmat([Pulse1Time; Pulse2Time; Pulse3Time], 1, length(RadiusIce));
        PulsePower = repmat([Pulse1Power; Pulse2Power; Pulse3Power], 1, length(RadiusIce));
        else
        PulseOnTime = repmat([Pulse1Time; Pulse2Time; Pulse3Time; Pulse4Time], 1, length(RadiusIce));
        PulsePower = repmat([Pulse1Power; Pulse2Power; Pulse3Power; Pulse4Power], 1, length(RadiusIce));
        end
    
    
%Mode2 - Laser Scan with 3 pulses, or 4 pulses changing 2nd and 3rd
elseif (length(Pulse2Time)>1 || (length(Pulse3Time)>1 && length(Pulse4Time)==1)) && RelativeTimingsMode==0 && length(SecondLaserOnTime)<2
        %Generate Radius Values
        disp("Varying pulse times (absolute, not pulse 4)")
        [Pulse2Time, Pulse3Time] = AbsolutePulseTiming(Pulse2Time, Pulse3Time) ;  
        %Mesh
             [NumLayers, MaxMassDiff, Method, Regions, ...
                 RadiusCHEndZone, LayersVapour, LayersIce, LayersCH, LayersCH1, LayersCHEndZone,...
                 RatioVapour, RatioIce, RatioCH, RatioCH1, RatioCHEndZone, MassDiff, ZoneMass]...
                 = MeshSolver(RadiusVapour, RadiusIce, RadiusCH, DensityVapour, FoamOnly);
        %For each laser combo
        for i=1:length(Pulse2Time)
        %Write File
             FileWriter(RadiusVapour, RadiusIce, RadiusCH, RadiusCHEndZone, LayersVapour, LayersIce, LayersCH, LayersCH1, LayersCHEndZone, RatioVapour, RatioIce, RatioCH, RatioCH1, RatioCHEndZone, DensityVapour, i, Regions, FoamOnly,TNBurn,...
                  Pulse1Time, Pulse2Time(i), Pulse3Time(i), Pulse4Time, Pulse1Power,Pulse2Power,Pulse3Power, Pulse4Power, RiseTime, LaserEndTime, ChangeTimeAt, TimeSpacing, ElectronHeatingOn, HeatingTimings, Heating, HeatingZones, HeatingEnergy, ...
     SecondLaser, SecondLaserOnTime, SecondLaserOffTime, SecondLaserWavelength, SecondLaserPower, Wavelength);
        end
    figure
    plot(MassDiff)
    ylabel('Percentage Mass Difference')
    yyaxis right;
    plot(ZoneMass, ':');
    ylabel('Zone Mass')
    title('Zoning Plot')
    xlabel('Zone')
    figure
    plot(Pulse3Time)
    hold on
    plot(Pulse2Time)
    hold off
    if isnan(Pulse4Time)==1
    PulseOnTime = [repmat(Pulse1Time, 1, length(Pulse2Time)); Pulse2Time; Pulse3Time];
    PulsePower = repmat([Pulse1Power; Pulse2Power; Pulse3Power], 1, length(Pulse2Time));
    else
    PulseOnTime = [repmat(Pulse1Time, 1, length(Pulse2Time)); Pulse2Time; Pulse3Time; repmat(Pulse4Time, 1, length(Pulse2Time))];
    PulsePower = repmat([Pulse1Power; Pulse2Power; Pulse3Power; Pulse4Power], 1, length(Pulse2Time));
    end
    
    
    %Mode2 - Laser Scan changing 3rd and 4th
elseif length(Pulse4Time)>1  && RelativeTimingsMode==0 && length(SecondLaserOnTime)<2
     disp("Varying pulse times (absolute, pulse 3 and 4)")
        %Generate Radius Values
        [Pulse3Time, Pulse4Time] = AbsolutePulseTiming(Pulse3Time, Pulse4Time) ;  
        %Mesh
             [NumLayers, MaxMassDiff, Method, Regions, ...
                 RadiusCHEndZone, LayersVapour, LayersIce, LayersCH, LayersCH1, LayersCHEndZone,...
                 RatioVapour, RatioIce, RatioCH, RatioCH1, RatioCHEndZone, MassDiff, ZoneMass]...
                 = MeshSolver(RadiusVapour, RadiusIce, RadiusCH, DensityVapour, FoamOnly);
        %For each laser combo
        for i=1:length(Pulse3Time)
        %Write File
             FileWriter(RadiusVapour, RadiusIce, RadiusCH, RadiusCHEndZone, LayersVapour, LayersIce, LayersCH, LayersCH1, LayersCHEndZone, RatioVapour, RatioIce, RatioCH, RatioCH1, RatioCHEndZone, DensityVapour, i, Regions, FoamOnly,TNBurn,...
                  Pulse1Time, Pulse2Time, Pulse3Time(i), Pulse4Time(i), Pulse1Power,Pulse2Power,Pulse3Power, Pulse4Power, RiseTime, LaserEndTime, ChangeTimeAt, TimeSpacing, ElectronHeatingOn, HeatingTimings, Heating, HeatingZones, HeatingEnergy, ...
     SecondLaser, SecondLaserOnTime, SecondLaserOffTime, SecondLaserWavelength, SecondLaserPower, Wavelength);
        end
    figure
    plot(MassDiff)
    ylabel('Percentage Mass Difference')
    yyaxis right;
    plot(ZoneMass, ':');
    ylabel('Zone Mass')
    title('Zoning Plot')
    xlabel('Zone')
    figure
    plot(Pulse4Time)
    hold on
    plot(Pulse3Time)
    hold off
    if isnan(Pulse4Time)==1
    PulseOnTime = [repmat(Pulse1Time, 1, length(Pulse3Time)); repmat(Pulse2Time, 1, length(Pulse4Time)); Pulse3Time; Pulse4Time];
    PulsePower = repmat([Pulse1Power; Pulse2Power; Pulse3Power; Pulse4Power], 1, length(Pulse4Time));
    else
    PulseOnTime = [repmat(Pulse1Time, 1, length(Pulse4Time)); repmat(Pulse2Time, 1, length(Pulse4Time)); Pulse3Time; Pulse4Time];
    PulsePower = repmat([Pulse1Power; Pulse2Power; Pulse3Power; Pulse4Power], 1, length(Pulse4Time));
    end
    
%Mode2.5 - Laser Scan with 4 pulses, changing 3rd and 4th pulses
elseif RelativeTimingsMode==1 && length(SecondLaserOnTime)<2
    disp("Varying pulse times (relative)")
    if length(Pulse4Time)==1
        Pulse4Time = repmat(Pulse4Time, 1, length(Pulse3Time));
    end
    if length(Pulse3Time)==1
        Pulse3Time = repmat(Pulse3Time, 1, length(Pulse4Time));
    end
    if length(Pulse2Time)==1
        Pulse2Time = repmat(Pulse2Time, 1, length(Pulse3Time));
    end
        %Mesh
             [NumLayers, MaxMassDiff, Method, Regions, ...
                 RadiusCHEndZone, LayersVapour, LayersIce, LayersCH, LayersCH1, LayersCHEndZone,...
                 RatioVapour, RatioIce, RatioCH, RatioCH1, RatioCHEndZone, MassDiff, ZoneMass]...
                 = MeshSolver(RadiusVapour, RadiusIce, RadiusCH, DensityVapour, FoamOnly);
        %For each laser combo
        for i=1:length(Pulse4Time)
        %Write File
             FileWriter(RadiusVapour, RadiusIce, RadiusCH, RadiusCHEndZone, LayersVapour, LayersIce, LayersCH, LayersCH1, LayersCHEndZone, RatioVapour, RatioIce, RatioCH, RatioCH1, RatioCHEndZone, DensityVapour, i, Regions, FoamOnly,TNBurn,...
                  Pulse1Time, Pulse2Time(i), Pulse3Time(i), Pulse4Time(i), Pulse1Power,Pulse2Power,Pulse3Power, Pulse4Power, RiseTime, LaserEndTime, ChangeTimeAt, TimeSpacing, ElectronHeatingOn, HeatingTimings, Heating, HeatingZones, HeatingEnergy, ...
     SecondLaser, SecondLaserOnTime, SecondLaserOffTime, SecondLaserWavelength, SecondLaserPower, Wavelength);
        end
    figure
    plot(MassDiff)
    ylabel('Percentage Mass Difference')
    yyaxis right;
    plot(ZoneMass, ':');
    ylabel('Zone Mass')
    title('Zoning Plot')
    xlabel('Zone')
    figure
    plot(Pulse3Time)
    hold on
    plot(Pulse2Time)
    plot(Pulse4Time)
    hold off
    if isnan(Pulse4Time(1))==1
    PulseOnTime = [repmat(Pulse1Time, 1, length(Pulse2Time)); Pulse2Time; Pulse3Time];
    PulsePower = repmat([Pulse1Power; Pulse2Power; Pulse3Power], 1, length(Pulse2Time));
    else
    PulseOnTime = [repmat(Pulse1Time, 1, length(Pulse2Time)); Pulse2Time; Pulse3Time; Pulse4Time];
    PulsePower = repmat([Pulse1Power; Pulse2Power; Pulse3Power; Pulse4Power], 1, length(Pulse2Time));
    end
    
    
    
%Mode=3 - Changing Electron Heating
elseif ElectronHeatingOn==1
    disp("Varying auxiliary heating")
    ArrayIndex=1;
    %Mesh
    [NumLayers, MaxMassDiff, Method, Regions, ...
         RadiusCHEndZone, LayersVapour, LayersIce, LayersCH, LayersCH1, LayersCHEndZone,...
         RatioVapour, RatioIce, RatioCH, RatioCH1, RatioCHEndZone, MassDiff, ZoneMass]...
         = MeshSolver(RadiusVapour, RadiusIce, RadiusCH, DensityVapour, FoamOnly);
     for i=1:size(HeatingTimings, 1)
        %Write File
             FileWriter(RadiusVapour, RadiusIce, RadiusCH, RadiusCHEndZone, LayersVapour, LayersIce, LayersCH, LayersCH1, LayersCHEndZone, RatioVapour, RatioIce, RatioCH, RatioCH1, RatioCHEndZone, DensityVapour, i, Regions, FoamOnly,TNBurn,...
                  Pulse1Time, Pulse2Time, Pulse3Time, Pulse4Time, Pulse1Power,Pulse2Power,Pulse3Power, Pulse4Power, RiseTime, LaserEndTime, ChangeTimeAt(i,:), TimeSpacing(i,:), ElectronHeatingOn, HeatingTimings(i,:), Heating(i,:), HeatingZones(i), HeatingEnergy(i),...
     SecondLaser, SecondLaserOnTime, SecondLaserOffTime, SecondLaserWavelength, SecondLaserPower, Wavelength);
        end
    %Plot Meshing
    figure
    plot(MassDiff)
    ylabel('Percentage Mass Difference')
    yyaxis right;
    plot(ZoneMass, ':');
    ylabel('Zone Mass')
    title('Zoning Plot')
    xlabel('Zone')
    %PulseOnTime and PulsePower
    if isnan(Pulse4Time)==1
    PulseOnTime = repmat([Pulse1Time; Pulse2Time; Pulse3Time], 1, size(HeatingTimings, 1));
    PulsePower = repmat([Pulse1Power; Pulse2Power; Pulse3Power], 1, size(HeatingTimings, 1));
    else
    PulseOnTime = repmat([Pulse1Time; Pulse2Time; Pulse3Time; Pulse4Time], 1, size(HeatingTimings, 1));
    PulsePower = repmat([Pulse1Power; Pulse2Power; Pulse3Power; Pulse4Power], 1, size(HeatingTimings, 1));
    end
    
%Mode=4 - ChangingLaserEndTime
elseif length(LaserEndTime)>1 
    disp("Varying laser end time only")
    ArrayIndex=1;
    %Mesh
    [NumLayers, MaxMassDiff, Method, Regions, ...
         RadiusCHEndZone, LayersVapour, LayersIce, LayersCH, LayersCH1, LayersCHEndZone,...
         RatioVapour, RatioIce, RatioCH, RatioCH1, RatioCHEndZone, MassDiff, ZoneMass]...
         = MeshSolver(RadiusVapour, RadiusIce, RadiusCH, DensityVapour, FoamOnly);
     for i=1:length(LaserEndTime)
        %Write File
             FileWriter(RadiusVapour, RadiusIce, RadiusCH, RadiusCHEndZone, LayersVapour, LayersIce, LayersCH, LayersCH1, LayersCHEndZone, RatioVapour, RatioIce, RatioCH, RatioCH1, RatioCHEndZone, DensityVapour, i, Regions, FoamOnly,TNBurn,...
                  Pulse1Time, Pulse2Time, Pulse3Time, Pulse4Time, Pulse1Power,Pulse2Power,Pulse3Power, Pulse4Power, RiseTime, LaserEndTime(i), ChangeTimeAt, TimeSpacing, ElectronHeatingOn, HeatingTimings, Heating, HeatingZones, HeatingEnergy, ...
     SecondLaser, SecondLaserOnTime, LaserEndTime(i), SecondLaserWavelength, SecondLaserPower, Wavelength);
        end
    %Plot Meshing
    figure
    plot(MassDiff)
    ylabel('Percentage Mass Difference')
    yyaxis right;
    plot(ZoneMass, ':');
    ylabel('Zone Mass')
    title('Zoning Plot')
    xlabel('Zone')
    %PulseOnTime and PulsePower
    if isnan(Pulse4Time)==1
    PulseOnTime = repmat([Pulse1Time; Pulse2Time; Pulse3Time], 1, length(LaserEndTime));
    PulsePower = repmat([Pulse1Power; Pulse2Power; Pulse3Power], 1, length(LaserEndTime));
    else
    PulseOnTime = repmat([Pulse1Time; Pulse2Time; Pulse3Time; Pulse4Time], 1, length(LaserEndTime));
    PulsePower = repmat([Pulse1Power; Pulse2Power; Pulse3Power; Pulse4Power], 1, length(LaserEndTime));
    end
    
    
    %Mode=5 - ChangingVapourPressure
elseif length(DensityVapour)>1
        disp('Changing vapour density')
    ArrayIndex=1;
     for i=1:length(DensityVapour)
         try
            %Mesh
             [NumLayers(i), MaxMassDiff(i), Method, Regions, ...
                 RadiusCHEndZone, LayersVapour, LayersIce, LayersCH, LayersCH1, LayersCHEndZone,...
                 RatioVapour, RatioIce, RatioCH, RatioCH1, RatioCHEndZone, ~, ~]...
                 = MeshSolver(RadiusVapour, RadiusIce, RadiusCH, DensityVapour(i), FoamOnly);
             %Write File
             FileWriter(RadiusVapour, RadiusIce, RadiusCH, RadiusCHEndZone, LayersVapour, LayersIce, LayersCH, LayersCH1, LayersCHEndZone, RatioVapour, RatioIce, RatioCH, RatioCH1, RatioCHEndZone, DensityVapour(i), i, Regions, FoamOnly,TNBurn,...
                  Pulse1Time, Pulse2Time, Pulse3Time, Pulse4Time, Pulse1Power,Pulse2Power,Pulse3Power, Pulse4Power, RiseTime, LaserEndTime, ChangeTimeAt, TimeSpacing, ElectronHeatingOn, HeatingTimings, Heating, HeatingZones, HeatingEnergy, ...
     SecondLaser, SecondLaserOnTime, SecondLaserOffTime, SecondLaserWavelength, SecondLaserPower, Wavelength);
              disp(['Meshing for file ', num2str(i), ' complete, number of layers = ', num2str(NumLayers(i)), ', max mass diff = ', num2str(MaxMassDiff(i))]);
         catch 
             warning('Meshing Error');
         end
        end
        figure
         plot(DensityVapour)
        figure
        yyaxis left
        plot(NumLayers)
        yyaxis right
        plot(MaxMassDiff)
        if isnan(Pulse4Time)==1
        PulseOnTime = repmat([Pulse1Time; Pulse2Time; Pulse3Time], 1, length(DensityVapour));
        PulsePower = repmat([Pulse1Power; Pulse2Power; Pulse3Power], 1, length(DensityVapour));
        else
        PulseOnTime = repmat([Pulse1Time; Pulse2Time; Pulse3Time; Pulse4Time], 1, length(DensityVapour));
        PulsePower = repmat([Pulse1Power; Pulse2Power; Pulse3Power; Pulse4Power], 1, length(DensityVapour));
        end
    
    
        %Mode=4 - Should work if only second laser on time is being changed
elseif RelativeTimingsMode==0  && length(SecondLaserOnTime)>1 && length(LaserEndTime)==1 && length(Pulse2Time)==1 && length(Pulse3Time)==1 && length(Pulse4Time)==1
    disp('Changing second laser time only')
    ArrayIndex=1;
    %Mesh
    [NumLayers, MaxMassDiff, Method, Regions, ...
         RadiusCHEndZone, LayersVapour, LayersIce, LayersCH, LayersCH1, LayersCHEndZone,...
         RatioVapour, RatioIce, RatioCH, RatioCH1, RatioCHEndZone, MassDiff, ZoneMass]...
         = MeshSolver(RadiusVapour, RadiusIce, RadiusCH, DensityVapour, FoamOnly);
     for i=1:length(SecondLaserOnTime)
        %Write File
             FileWriter(RadiusVapour, RadiusIce, RadiusCH, RadiusCHEndZone, LayersVapour, LayersIce, LayersCH, LayersCH1, LayersCHEndZone, RatioVapour, RatioIce, RatioCH, RatioCH1, RatioCHEndZone, DensityVapour, i, Regions, FoamOnly,TNBurn,...
                  Pulse1Time, Pulse2Time, Pulse3Time, Pulse4Time, Pulse1Power,Pulse2Power,Pulse3Power, Pulse4Power, RiseTime, LaserEndTime, ChangeTimeAt, TimeSpacing, ElectronHeatingOn, HeatingTimings, Heating, HeatingZones, HeatingEnergy, ...
     SecondLaser, SecondLaserOnTime(i), SecondLaserOffTime, SecondLaserWavelength, SecondLaserPower, Wavelength);
        end
    %Plot Meshing
    figure
    plot(MassDiff)
    ylabel('Percentage Mass Difference')
    yyaxis right;
    plot(ZoneMass, ':');
    ylabel('Zone Mass')
    title('Zoning Plot')
    xlabel('Zone')
    %PulseOnTime and PulsePower
    if isnan(Pulse4Time)==1
    PulseOnTime = repmat([Pulse1Time; Pulse2Time; Pulse3Time], 1, length(SecondLaserOnTime));
    PulsePower = repmat([Pulse1Power; Pulse2Power; Pulse3Power], 1, length(SecondLaserOnTime));
    else
    PulseOnTime = repmat([Pulse1Time; Pulse2Time; Pulse3Time; Pulse4Time], 1, length(SecondLaserOnTime));
    PulsePower = repmat([Pulse1Power; Pulse2Power; Pulse3Power; Pulse4Power], 1, length(SecondLaserOnTime));
    end
    
    
    
   
    
    %Changing second laser on time and pulse 4 time (absolute)
elseif length(Pulse4Time)>1  && RelativeTimingsMode==0 && length(SecondLaserOnTime)>1 && length(LaserEndTime)==1
        disp('Changing second laser time and pulse 4 (absolute)')
        %Generate Radius Values
        [Pulse4Time, SecondLaserOnTime] = AbsolutePulseTiming(Pulse4Time, SecondLaserOnTime) ;  
        %Mesh
             [NumLayers, MaxMassDiff, Method, Regions, ...
                 RadiusCHEndZone, LayersVapour, LayersIce, LayersCH, LayersCH1, LayersCHEndZone,...
                 RatioVapour, RatioIce, RatioCH, RatioCH1, RatioCHEndZone, MassDiff, ZoneMass]...
                 = MeshSolver(RadiusVapour, RadiusIce, RadiusCH, DensityVapour, FoamOnly);
        %For each laser combo
        for i=1:length(Pulse4Time)
        %Write File
             FileWriter(RadiusVapour, RadiusIce, RadiusCH, RadiusCHEndZone, LayersVapour, LayersIce, LayersCH, LayersCH1, LayersCHEndZone, RatioVapour, RatioIce, RatioCH, RatioCH1, RatioCHEndZone, DensityVapour, i, Regions, FoamOnly,TNBurn,...
                  Pulse1Time, Pulse2Time, Pulse3Time, Pulse4Time(i), Pulse1Power,Pulse2Power,Pulse3Power, Pulse4Power, RiseTime, LaserEndTime, ChangeTimeAt, TimeSpacing, ElectronHeatingOn, HeatingTimings, Heating, HeatingZones, HeatingEnergy, ...
     SecondLaser, SecondLaserOnTime(i), SecondLaserOffTime, SecondLaserWavelength, SecondLaserPower, Wavelength);
        end
    figure
    plot(MassDiff)
    ylabel('Percentage Mass Difference')
    yyaxis right;
    plot(ZoneMass, ':');
    ylabel('Zone Mass')
    title('Zoning Plot')
    xlabel('Zone')
    figure
    plot(Pulse4Time)
    hold on
    plot(SecondLaserOnTime)
    hold off
    if isnan(Pulse4Time)==1
    PulseOnTime = [repmat(Pulse1Time, 1, length(Pulse3Time)); repmat(Pulse2Time, 1, length(Pulse4Time)); repmat(Pulse3Time, 1, length(Pulse4Time)); Pulse4Time];
    PulsePower = repmat([Pulse1Power; Pulse2Power; Pulse3Power; Pulse4Power], 1, length(Pulse4Time));
    else
    PulseOnTime = [repmat(Pulse1Time, 1, length(Pulse4Time)); repmat(Pulse2Time, 1, length(Pulse4Time)); repmat(Pulse3Time, 1, length(Pulse4Time)); Pulse4Time];
    PulsePower = repmat([Pulse1Power; Pulse2Power; Pulse3Power; Pulse4Power], 1, length(Pulse4Time));
    end
    
     %Changing second laser on time and pulse 3 time (absolute)
    elseif length(Pulse3Time)>1  && RelativeTimingsMode==0 && length(SecondLaserOnTime)>1 && length(LaserEndTime)==1
        disp('Changing second laser time and pulse 3 (absolute)')
        %Generate Radius Values
        [Pulse3Time, SecondLaserOnTime] = AbsolutePulseTiming(Pulse3Time, SecondLaserOnTime) ;  
        %Mesh
             [NumLayers, MaxMassDiff, Method, Regions, ...
                 RadiusCHEndZone, LayersVapour, LayersIce, LayersCH, LayersCH1, LayersCHEndZone,...
                 RatioVapour, RatioIce, RatioCH, RatioCH1, RatioCHEndZone, MassDiff, ZoneMass]...
                 = MeshSolver(RadiusVapour, RadiusIce, RadiusCH, DensityVapour, FoamOnly);
        %For each laser combo
        for i=1:length(Pulse3Time)
        %Write File
             FileWriter(RadiusVapour, RadiusIce, RadiusCH, RadiusCHEndZone, LayersVapour, LayersIce, LayersCH, LayersCH1, LayersCHEndZone, RatioVapour, RatioIce, RatioCH, RatioCH1, RatioCHEndZone, DensityVapour, i, Regions, FoamOnly,TNBurn,...
                  Pulse1Time, Pulse2Time, Pulse3Time(i), Pulse4Time, Pulse1Power,Pulse2Power,Pulse3Power, Pulse4Power, RiseTime, LaserEndTime, ChangeTimeAt, TimeSpacing, ElectronHeatingOn, HeatingTimings, Heating, HeatingZones, HeatingEnergy, ...
     SecondLaser, SecondLaserOnTime(i), SecondLaserOffTime, SecondLaserWavelength, SecondLaserPower, Wavelength);
        end
    figure
    plot(MassDiff)
    ylabel('Percentage Mass Difference')
    yyaxis right;
    plot(ZoneMass, ':');
    ylabel('Zone Mass')
    title('Zoning Plot')
    xlabel('Zone')
    figure
    plot(Pulse3Time)
    hold on
    plot(SecondLaserOnTime)
    hold off
    if isnan(Pulse4Time)==1
    PulseOnTime = [repmat(Pulse1Time, 1, length(Pulse3Time)); repmat(Pulse2Time, 1, length(Pulse3Time)); Pulse3Time; repmat(Pulse4Time, 1, length(Pulse3Time))];
    PulsePower = repmat([Pulse1Power; Pulse2Power; Pulse3Power; Pulse4Power], 1, length(Pulse3Time));
    else
    PulseOnTime = [repmat(Pulse1Time, 1, length(Pulse3Time)); repmat(Pulse2Time, 1, length(Pulse3Time)); Pulse3Time; repmat(Pulse4Time, 1, length(Pulse3Time))];
    PulsePower = repmat([Pulse1Power; Pulse2Power; Pulse3Power; Pulse4Power], 1, length(Pulse3Time));
    end
    
     %Changing second laser on time and pulse 2 time (absolute)
    elseif length(Pulse2Time)>1  && RelativeTimingsMode==0 && length(SecondLaserOnTime)>1 && length(LaserEndTime)==1 
        disp('Changing second laser time and pulse 2 (absolute)')
        %Generate Radius Values
        [Pulse2Time, SecondLaserOnTime] = AbsolutePulseTiming(Pulse2Time, SecondLaserOnTime) ;  
        %Mesh
             [NumLayers, MaxMassDiff, Method, Regions, ...
                 RadiusCHEndZone, LayersVapour, LayersIce, LayersCH, LayersCH1, LayersCHEndZone,...
                 RatioVapour, RatioIce, RatioCH, RatioCH1, RatioCHEndZone, MassDiff, ZoneMass]...
                 = MeshSolver(RadiusVapour, RadiusIce, RadiusCH, DensityVapour, FoamOnly);
        %For each laser combo
        for i=1:length(Pulse2Time)
        %Write File
             FileWriter(RadiusVapour, RadiusIce, RadiusCH, RadiusCHEndZone, LayersVapour, LayersIce, LayersCH, LayersCH1, LayersCHEndZone, RatioVapour, RatioIce, RatioCH, RatioCH1, RatioCHEndZone, DensityVapour, i, Regions, FoamOnly,TNBurn,...
                  Pulse1Time, Pulse2Time(i), Pulse3Time, Pulse4Time, Pulse1Power,Pulse2Power,Pulse3Power, Pulse4Power, RiseTime, LaserEndTime, ChangeTimeAt, TimeSpacing, ElectronHeatingOn, HeatingTimings, Heating, HeatingZones, HeatingEnergy, ...
     SecondLaser, SecondLaserOnTime(i), SecondLaserOffTime, SecondLaserWavelength, SecondLaserPower, Wavelength);
        end
    figure
    plot(MassDiff)
    ylabel('Percentage Mass Difference')
    yyaxis right;
    plot(ZoneMass, ':');
    ylabel('Zone Mass')
    title('Zoning Plot')
    xlabel('Zone')
    figure
    plot(Pulse2Time)
    hold on
    plot(SecondLaserOnTime)
    hold off
    if isnan(Pulse4Time)==1
    PulseOnTime = [repmat(Pulse1Time, 1, length(Pulse2Time)); Pulse2Time; repmat(Pulse3Time, 1, length(Pulse2Time)); repmat(Pulse4Time, 1, length(Pulse2Time))];
    PulsePower = repmat([Pulse1Power; Pulse2Power; Pulse3Power; Pulse4Power], 1, length(Pulse2Time));
    else
    PulseOnTime = [repmat(Pulse1Time, 1, length(Pulse2Time)); Pulse2Time; repmat(Pulse3Time, 1, length(Pulse2Time)); repmat(Pulse4Time, 1, length(Pulse2Time))];
    PulsePower = repmat([Pulse1Power; Pulse2Power; Pulse3Power; Pulse4Power], 1, length(Pulse2Time));
    end
    
    
 %Changing second laser on time any other pulse (relative)
elseif RelativeTimingsMode==1  && length(SecondLaserOnTime)>1
    disp('Changing second laser time and other pulse times (relative)')
    if length(Pulse4Time)==1 
        Pulse4Time = repmat(Pulse4Time, 1, length(SecondLaserOnTime));
    end
    if length(Pulse3Time)==1 
        Pulse3Time = repmat(Pulse3Time, 1, length(SecondLaserOnTime));
    end
    if length(Pulse2Time)==1 
        Pulse2Time = repmat(Pulse2Time, 1, length(SecondLaserOnTime));
    end
            %Mesh
             [NumLayers, MaxMassDiff, Method, Regions, ...
                 RadiusCHEndZone, LayersVapour, LayersIce, LayersCH, LayersCH1, LayersCHEndZone,...
                 RatioVapour, RatioIce, RatioCH, RatioCH1, RatioCHEndZone, MassDiff, ZoneMass]...
                 = MeshSolver(RadiusVapour, RadiusIce, RadiusCH, DensityVapour, FoamOnly);
        %For each laser combo
        for i=1:length(Pulse4Time)
        %Write File
             FileWriter(RadiusVapour, RadiusIce, RadiusCH, RadiusCHEndZone, LayersVapour, LayersIce, LayersCH, LayersCH1, LayersCHEndZone, RatioVapour, RatioIce, RatioCH, RatioCH1, RatioCHEndZone, DensityVapour, i, Regions, FoamOnly,TNBurn,...
                  Pulse1Time, Pulse2Time(i), Pulse3Time(i), Pulse4Time(i), Pulse1Power,Pulse2Power,Pulse3Power, Pulse4Power, RiseTime, LaserEndTime, ChangeTimeAt, TimeSpacing, ElectronHeatingOn, HeatingTimings, Heating, HeatingZones, HeatingEnergy, ...
     SecondLaser, SecondLaserOnTime(i), SecondLaserOffTime, SecondLaserWavelength, SecondLaserPower, Wavelength);
        end
    figure
    plot(MassDiff)
    ylabel('Percentage Mass Difference')
    yyaxis right;
    plot(ZoneMass, ':');
    ylabel('Zone Mass')
    title('Zoning Plot')
    xlabel('Zone')
    figure
    plot(Pulse3Time)
    hold on
    plot(Pulse2Time)
    plot(Pulse4Time)
    plot(SecondLaserOnTime)
    hold off
    if isnan(Pulse4Time(1))==1
    PulseOnTime = [repmat(Pulse1Time, 1, length(Pulse2Time)); Pulse2Time; Pulse3Time];
    PulsePower = repmat([Pulse1Power; Pulse2Power; Pulse3Power], 1, length(Pulse2Time));
    else
    PulseOnTime = [repmat(Pulse1Time, 1, length(Pulse2Time)); Pulse2Time; Pulse3Time; Pulse4Time];
    PulsePower = repmat([Pulse1Power; Pulse2Power; Pulse3Power; Pulse4Power], 1, length(Pulse2Time));
    end
    
      
        
%Mode=0 - Single Run
    else
    disp('Single run only')
    ArrayIndex=1;
    %Mesh
    [NumLayers, MaxMassDiff, Method, Regions, ...
         RadiusCHEndZone, LayersVapour, LayersIce, LayersCH, LayersCH1, LayersCHEndZone,...
         RatioVapour, RatioIce, RatioCH, RatioCH1, RatioCHEndZone, MassDiff, ZoneMass]...
         = MeshSolver(RadiusVapour, RadiusIce, RadiusCH, DensityVapour, FoamOnly);
    %Write File
    FileWriter(RadiusVapour, RadiusIce, RadiusCH, RadiusCHEndZone, LayersVapour, LayersIce, LayersCH, LayersCH1, LayersCHEndZone, RatioVapour, RatioIce, RatioCH, RatioCH1, RatioCHEndZone, DensityVapour, ArrayIndex, Regions,FoamOnly,TNBurn,...
         Pulse1Time, Pulse2Time, Pulse3Time, Pulse4Time, Pulse1Power,Pulse2Power,Pulse3Power, Pulse4Power, RiseTime, LaserEndTime, ChangeTimeAt, TimeSpacing, ElectronHeatingOn, HeatingTimings, Heating, HeatingZones, HeatingEnergy, ...
     SecondLaser, SecondLaserOnTime, SecondLaserOffTime, SecondLaserWavelength, SecondLaserPower, Wavelength);
    %Plot Meshing
    figure
    plot(MassDiff)
    ylabel('Percentage Mass Difference')
    yyaxis right;
    plot(ZoneMass, ':');
    ylabel('Zone Mass')
    title('Zoning Plot')
    xlabel('Zone')
    %PulseOnTime and PulsePower
    PulseOnTime = repmat([Pulse1Time; Pulse2Time; Pulse3Time], 1, length(RadiusIce));
    PulsePower = repmat([Pulse1Power; Pulse2Power; Pulse3Power], 1, length(RadiusIce));
    end
    
     %Calculate time length    
    Time=[0];
    for i=1:(length(ChangeTimeAt)-1)
        Time = [Time, (ChangeTimeAt(i)+TimeSpacing(i)):TimeSpacing(i):ChangeTimeAt(i+1)];
    end
        Timelength = length(Time)

 %Save the change variables to a file, and back up this script in save
    %location (so if inf files are deleted for space it's easy to reproduce
    %them)
    if ElectronHeatingOn==1
       save('MultiFile/Variables','LaserEndTime', 'PulsePower', 'PulseOnTime', 'HeatingTimings', 'Heating', 'HeatingZones', 'HeatingEnergy')
    elseif SecondLaser==1
        save('MultiFile/Variables','LaserEndTime', 'PulsePower', 'PulseOnTime', 'SecondLaserOnTime', 'SecondLaserOffTime', 'SecondLaserWavelength', 'SecondLaserPower')
    else
    save('MultiFile/Variables','LaserEndTime', 'PulsePower', 'PulseOnTime')
    end
    FileNameAndLocation=[mfilename('fullpath')];
    newbackup=sprintf('%sGeneratorbackup.m','MultiFile/');
    currentfile=strcat(FileNameAndLocation, '.m');
    copyfile(currentfile,newbackup);
        
end

function [NumLayers, MaxMassDiff, Method, Regions, ...
    RadiusCHEndZone, LayersVapour, LayersIce, LayersCH, LayersCH1, LayersCHEndZone,...
    RatioVapour, RatioIce, RatioCH, RatioCH1, RatioCHEndZone, MassDiff, ZoneMass]...
    = MeshSolver(RadiusVapour, RadiusIce, RadiusCH, DensityVapour, FoamOnly)
Meshing_Git
MaxMassDiff = max(abs(MassDiff(LayersVapour:end)));
end

function [] = FileWriter(RadiusVapour, RadiusIce, RadiusCH, RadiusCHEndZone, LayersVapour, LayersIce, LayersCH, LayersCH1, LayersCHEndZone, RatioVapour, RatioIce, RatioCH, RatioCH1, RatioCHEndZone, DensityVapour, ArrayIndex, Regions, FoamOnly, TNBurn, ...
    Pulse1Time, Pulse2Time, Pulse3Time, Pulse4Time, Pulse1Power,Pulse2Power,Pulse3Power, Pulse4Power, RiseTime, LaserEndTime, ChangeTimeAt, TimeSpacing, ElectronHeatingOn, HeatingTimings, Heating, HeatingZones, HeatingEnergy, ...
     SecondLaser, SecondLaserOnTime, SecondLaserOffTime, SecondLaserWavelength, SecondLaserPower, Wavelength)
if isnan(Pulse4Time)
    clear Pulse4Time
end
BatchFileWriterAdaptive_Git
end
