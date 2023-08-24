%This script allows the user to specify a range of different settings for
%an ICF implosion, which are then used to create input decks for Hyades. It
%has a range of functionality, and can be used for batch runs.

%A parameter "Selection" allows the user to have multiple sets of settings
%saved at once, and choose between them.

%Implosions can be done at any wavelength. They can be 'two-colour'
%implosions at another wavelength, and/or include auxiliary heating. Foam
%only hydrodynamic equivalent capsules can also be used.

%Regardless of these settings, batch lots of input decks can be produced
%whether different parameters are varied. These options are outlined below:

%Change thicknesses - the desired vapour and ice boundaries are saved as
%arrays. Each permutation of these two arrays will be produced (each
%capsule must be meshed individually).

%Change laser timings (RelativeTimingMode OFF) - set two different laser
%timings in the arrays. All possible combinations will be ran. Arrays can
%be set as constants relative to each other (i.e. pulse 4 = pulse 3 + x).

%Change laser timings (RelativeTimingMode ON) - set the different laser
%timings in the arrays. No combinations will be performed: the arrays must
%be the same length, and a run will be performed using all the 1st values,
%then all the second values, etc. All the permutations of two arrays can be
%calculated using AbsolutePulseTimings. One array can be varied relative to
%another using RelativePulseTimings (and the first will be changed so that
%the arrays have the same lengths and pairings). If one of the arrays is a
%single constant value, this will be automatically replciated to the other
%array lengths.

%The second laser can be turned on using Second Laser parameter.
%SecondLaserOnTime is treated just like any other pulse time, and can be
%varied in combination with the others. SecondLaserOffTime is a legacy
%variable - this should just equal laser end time. 

%Auxiliary heating can also be applied. Here, you should first specify the
%file for the unheated capsule, and make sure the radius/laser inputs use
%match that input deck. The timing/magnitude of heating can then be varied.

%Other parameters (laser end times, vapour densities) can also be varied
%(seperately from other variables) by feeding them arrays as inputs.

%FoamOnly makes the second layer a dry CD foam (no wetting). TNBurn enables
%the thermonuclear burn in layers 1 and 2. Electron Heating and Second
%Laser turn these two things on.

%SizeMod is used to ensure the laser power scales correctly with size.
%Rather than changing radius directly, scale it to a 0.285cm capsule
%(sizemod=1) and multiply by sizemod. Sizemod is then used to adjust the
%laser power accordingly. 
 
%I have added a series of different examples (saved as 'selections' below
%to try and demonstrate the different functionality).

%The files are saved in the folder 'MultiFile'. You many need to create
%such a folder. They will be produced as File1 to FileN. If these files
%already exist, they will be overwritten. If File 50 already exists and you
%produce 30 new files, File1 to File30 will be overwritten - but File40 and
%File50 will be unchanged from the previous run.


%% Use Selection to choose between different examples - or create your own!
Selection=103;




%% Example 1 - Simple example producing a single input deck for a third
%harmonic implosion.
if Selection==101
    
    %Wavelength - Set wavelength of laser. 0.351 for two colour, 0.193 for
    %ArF
    Wavelength = 0.351;
    
    %SizeMod is required for every run. SizeMod = OuterRadius/2.85mm. This
    %is then used to scale e.g. laser powers, and for the meshing script.
     SizeMod = 0.75;
     
    %RelativeTimingsMode is used to tell the function how to vary certain
    %parameters. Should be 0, unless you are using the RelativeTimingsMode
    %functionality when varying laser timings.
    RelativeTimingsMode = 0;
    
    %FoamOnly=1 makes hydroequivalent capsules with no DT in foam layer.
    %Use FoamOnly=0 for normal implosions.
    FoamOnly=0;
    
    %TNBurn=1 activates the TN burn package (fusion reactions).
    TNBurn=1;
       
    %Set your pulse timings (i.e. switch on time of each pulse)
    Pulse1Time = 0;
    Pulse2Time = 3.6.*10^-9;
    Pulse3Time = 6.8.*10^-9;
    Pulse4Time = 8.2.*10^-9;

    %Pulse Powers are automatically set using sizemod and wavelength to
    %satisy our irradiance condition.
    Pulse1Power = 2*((0.351/Wavelength)^2)*SizeMod^2;
    Pulse2Power = 14*((0.351/Wavelength)^2)*SizeMod^2;;
    Pulse3Power = 98*((0.351/Wavelength)^2)*SizeMod^2;
    Pulse4Power = 692.13*((0.351/Wavelength)^2)*SizeMod^2;
    
    %Other laser settings. Rise time is how long any changes in laser power
    %take (0.2 ns throughout all my sims). Laser end time is when the laser
    %begins to be turned off.
    RiseTime = 2e-10;
    LaserEndTime =16.2.*10^-9;
    
    %Set the capsule dimensions (cm) and vapour density (g/cm^3)
    RadiusVapour = 0.19275;
    RadiusIce = 0.21;
    RadiusCH= 0.285 * SizeMod;
    DensityVapour= 0.0009;
    
    %Set time spacing for the simulation. At ChangeTime(i), the TimeSpacing
    %is set to TimeSpacing(i). The final time in ChangeTimeAt is the end
    %time for the simulation. For this reason, ChangeTimeAt should be one
    %element longer than TimeSpacing.
    ChangeTimeAt = [0 15 24].*10^-9;
    TimeSpacing = [1E-10 1E-11];
    
    %ElectronHeatingOn=0 for no electron heating - all other values =NaN
    %unless ElectronHeatingOn=1.
    ElectronHeatingOn=0;
    HeatingTimings = NaN;
    Heating = NaN;
    HeatingZones = NaN;
    HeatingEnergy = NaN;
        
    %SecondLaser=1 only if we are using the second laser. Like for electron
    %heating, other values are NaN unless SecondLaser=1.
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


%% Example 2 - Include the second laser, and produce a batch run of input
%decks where the thickness is varied.
if Selection==102
    
    %General settings as defined above   
    Wavelength = 0.351;
    SizeMod = 0.5;
    RelativeTimingsMode = 0;
    FoamOnly=0;
    TNBurn=1;
      
    %SecondLaserOnTime treated like any other pulse time
    Pulse1Time = 0;
    Pulse2Time = [0.9].*10^-9;
    Pulse3Time = [3.8].*10^-9;
    Pulse4Time = [6.2].*10^-9;
    SecondLaserOnTime = [9.2]*10^-9;

    Pulse1Power = 2*((0.351/Wavelength)^2)*SizeMod^2;
    Pulse2Power = 14*((0.351/Wavelength)^2)*SizeMod^2;;
    Pulse3Power = 98*((0.351/Wavelength)^2)*SizeMod^2;
    Pulse4Power = 692.13*((0.351/Wavelength)^2)*SizeMod^2;
    RiseTime = 2e-10;
    LaserEndTime = [11].*10^-9;
    
    %There are two arrays for thickness. An input deck will be produced for
    %every combination of these two arrays. (It is possible to include the
    %second laser while any variable is being changed).
    RadiusVapour = [0.123:0.0005:0.127];
    RadiusIce = [0.135:0.0005:0.139];
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
    %OffTime, Wavelength, Radius and Power should be single valued
    %(LaserOnTime can be changed like any other pulse time - but not in
    %this example, as thickness is being changed instead). SecondLaserFn
    %determines the appropriate power, given the wavelength and radius.
    SecondLaser=1;
    SecondLaserWavelength = 0.193; %ArF frequency
    SecondLaserOffTime = LaserEndTime; %Legacy variable - should equal LaserEndTime. 
    SecondLaserRadius = RadiusCH; %Legacy variable - should equal RadiusCH. 
    [SecondLaserPower] = SecondLaserFn(SecondLaserWavelength, SecondLaserRadius);
  
    %Same function is used - always unchanged.
    [NumLayers, MassDiff, Method, Regions, ...
    Time, PulseOnTime, PulsePower] = GeneratorFunction(...
    RadiusVapour, RadiusIce, RadiusCH, DensityVapour, FoamOnly, TNBurn,...
    Pulse1Time, Pulse2Time, Pulse3Time, Pulse4Time, Pulse1Power, Pulse2Power, Pulse3Power, Pulse4Power, ...
    RelativeTimingsMode, RiseTime, LaserEndTime, ChangeTimeAt, TimeSpacing,...
    ElectronHeatingOn, HeatingTimings, Heating, HeatingZones, HeatingEnergy,...
    SecondLaser, SecondLaserOnTime, LaserEndTime, SecondLaserWavelength, SecondLaserPower, Wavelength);
end



%% Example 3 - Include auxiliary heating, and change the timing this is
%applied at.
if Selection==103
    
    Wavelength = 0.351;
    SizeMod = 0.5;
    RelativeTimingsMode = 0;
    FoamOnly=0;
    TNBurn=1;
    
    %We now specify the file for the unheated capsule. This is used to
    %identify the bang time, and the vapour mass - so that the timings and
    %amount of heating energy are appropriate.
    UnheatedFile = '\\aldaq1.physics.ox.ac.uk\Archer\Robert\Desktop Hyades Files\Dated\230106\TwoColourPt5\0pt5TwoColour.cdf';
    
    Pulse1Time = 0;
    Pulse2Time = [1.95].*10^-9;
    Pulse3Time = [4.35].*10^-9;
    Pulse4Time = [8.05].*10^-9;
    SecondLaserOnTime = [10]*10^-9;

    Pulse1Power = 2*((0.351/Wavelength)^2)*SizeMod^2;
    Pulse2Power = 14*((0.351/Wavelength)^2)*SizeMod^2;;
    Pulse3Power = 98*((0.351/Wavelength)^2)*SizeMod^2;
    Pulse4Power = 692.13*((0.351/Wavelength)^2)*SizeMod^2;
    RiseTime = 2e-10;
    LaserEndTime = [12.8].*10^-9;
    
    RadiusVapour = [0.123];
    RadiusIce = [0.136];
    RadiusCH=0.1425;
    DensityVapour=[0.00102];
    ChangeTimeAt = [0 12 14 16].*10^-9;
    TimeSpacing = [1E-10 1E-11 1E-10];
    
    %We turn the auxiliary heating on, and create the following empty
    %variables (filled in the next chunk of code). The amount of heating
    %is specified in HeatingkJ. UnheatedFileReader reads the file to obtain the
    %necessary variables.
    ElectronHeatingOn=1;
    HeatingTimings = [];
    Heating = [];
    HeatingZones = [];
    ChangeTimeAt = [];
    TimeSpacing = [];
    HeatingEnergy = [];
    [BangTime, Time, TotalVapourMass, IceBoundaryIndex] = UnheatedFileReader(UnheatedFile)
    HeatingkJ = 4;

    %In the for loop, we specify the range of times around the bang time
    %that we wish to apply our heating at. We then identify these times and
    %work out the appropriate power, and fill in our empty arrays from
    %above.
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
    %Correct the units.
    HeatingTimings=HeatingTimings*10^-9;
    ChangeTimeAt=ChangeTimeAt*10^-9;
    HeatingEnergy = HeatingEnergy*1000;
    
    %We can use both second laser and auxiliary heating if we want.
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

%% Example 4 - Changing power of auxiliary heating
if Selection==104
    
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

    Pulse1Power = 2*((0.351/Wavelength)^2)*SizeMod^2;
    Pulse2Power = 14*((0.351/Wavelength)^2)*SizeMod^2;
    Pulse3Power = 98*((0.351/Wavelength)^2)*SizeMod^2;
    Pulse4Power = 692.13*((0.351/Wavelength)^2)*SizeMod^2;
    RiseTime = 2e-10;
    LaserEndTime = [11].*10^-9;
    
    RadiusVapour = [0.125];
    RadiusIce = [0.137];
    RadiusCH=0.1425;
    DensityVapour=[0.0012];
    ChangeTimeAt = [0 5 8 11 15].*10^-9;
    TimeSpacing = [1E-10 1E-11 1E-10 1E-11];
    
    %Electron heating turned on as in last example, and unheated file read. This time
    %we specify our fixed HeatingStartTime (likely the optimal found from
    %running a timing scan).
    ElectronHeatingOn=1;
    HeatingTimings = [];
    Heating = [];
    HeatingZones = [];
    ChangeTimeAt = [];
    TimeSpacing = [];
    HeatingEnergy = [];
    [BangTime, Time, TotalVapourMass, IceBoundaryIndex] = UnheatedFileReader(UnheatedFile)
    HeatingStartTime = 13.2552;
   
    %We now specify our heating energy as an array in the for loop.
    for HeatingkJ = [2:2:60] 
    [~, TimeIndex] = min(abs(Time - HeatingStartTime*10^-9)); %Total Vapour mass is actually the same regardless of time, but done this for continuity with above.
    HeatingConv = HeatingkJ*(10^10)/(5E-13 * TotalVapourMass(TimeIndex)); %EnergyKJ * KJ2Erg / (Time*Mass) to give Erg/g/sec
    HeatingTimings = [HeatingTimings; HeatingStartTime HeatingStartTime+0.0001 HeatingStartTime+0.0005 HeatingStartTime+0.0006];
    Heating = [Heating; 0 HeatingConv HeatingConv 0];
    HeatingZones = [HeatingZones; IceBoundaryIndex];
    ChangeTimeAt = [ChangeTimeAt; 0 11 HeatingStartTime HeatingStartTime+0.002 HeatingStartTime+0.02 Time(end).*10^9];
    TimeSpacing = [TimeSpacing;  1E-10 1E-11 1E-13 1E-12 2E-11 ];
    HeatingEnergy = [HeatingEnergy; HeatingkJ];
    end
    %Fix the units
    HeatingTimings=HeatingTimings*10^-9;
    ChangeTimeAt=ChangeTimeAt*10^-9;
    HeatingEnergy = HeatingEnergy*1000;
    
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


%% Example 5 - Changing laser timings (absolute timings mode)
if Selection==105
    
    Wavelength = 0.351;
    SizeMod = 0.5;
    RelativeTimingsMode = 0;
    FoamOnly=0;
    TNBurn=1;
    
    %Here, we change our laser timings. We specify two of our pulse times
    %to be arrays - and we will produce an input deck for each permutation.
    %Here I have set Pulse2Time to be a single value for all runs, but we
    %could use the commented line to also set this to change with
    %Pulse3Time without requiring any further changes to the settings/code.
    Pulse1Time = 0;
    Pulse2Time = [1.2].*10-9;
    Pulse3Time = [3.6:0.1:4.6].*10^-9;
    Pulse4Time = [5.5:0.1:6.2].*10^-9;  
%     Pulse2Time = Pulse3Time - 2.4*10^-9;
   
    Pulse1Power = 2*((0.351/Wavelength)^2)*SizeMod^2;
    Pulse2Power = 14*((0.351/Wavelength)^2)*SizeMod^2;
    Pulse3Power = 98*((0.351/Wavelength)^2)*SizeMod^2;
    Pulse4Power = 692.13*((0.351/Wavelength)^2)*SizeMod^2;
    RiseTime = 2e-10;
    LaserEndTime = [11].*10^-9;
    
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
    
    %Second Laser could be on or off - not using for this run.
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
    SecondLaser, SecondLaserOnTime, LaserEndTime, SecondLaserWavelength, SecondLaserPower, Wavelength);
end


%% Example 6 - Changing laser timings with RelativeTimingsMode.
if Selection==106
    
    %Note RelativeTimingsMode=1, as we will use Relative timing function
    Wavelength = 0.351;
    SizeMod = 0.5;
    RelativeTimingsMode = 1;
    FoamOnly=0;
    TNBurn=1;
    
    %Here, we use RelativePulseTiming function. First, we specify an array
    %of values for pulse 3. We then set SecondLaserOnTime to vary between 4
    %and 5ns after each value of Pulse 3 Time (with 0.1ns spacing). We have
    %also said that Pulse 2 Time will be 2.4ns earlier than pulse 3 time,
    %regardless of the value. We need to set RelativeTimingsMode=1 for this
    %to work.
    Pulse1Time = 0;
    Pulse3Time = [3.6:0.1:4.6].*10^-9;
    [Pulse3Time, SecondLaserOnTime] = RelativePulseTiming(Pulse3Time, 4, 5, 0.1) ;  
    Pulse2Time = Pulse3Time - 2.4*10^-9;
    Pulse4Time = [5]*10^-9;

    Pulse1Power = 2*((0.351/Wavelength)^2)*SizeMod^2;
    Pulse2Power = 14*((0.351/Wavelength)^2)*SizeMod^2;;
    Pulse3Power = 98*((0.351/Wavelength)^2)*SizeMod^2;
    Pulse4Power = 692.13*((0.351/Wavelength)^2)*SizeMod^2;
    RiseTime = 2e-10;
    LaserEndTime = [11].*10^-9;
    
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


%% Example 7 - Changing vapour density
if Selection==107
    
    %Here I have changed the wavelength, to give us an ArF implosion
    Wavelength = 0.193;
    SizeMod = 0.5;
    RelativeTimingsMode = 0;
    FoamOnly=0;
    TNBurn=1;
    
    Pulse1Time = 0;
    Pulse2Time = [0.9].*10^-9;
    Pulse3Time = [3.8].*10^-9;
    Pulse4Time = [6.2].*10^-9;
    
    Pulse1Power = 2*((0.351/Wavelength)^2)*SizeMod^2;
    Pulse2Power = 14*((0.351/Wavelength)^2)*SizeMod^2;;
    Pulse3Power = 98*((0.351/Wavelength)^2)*SizeMod^2;
    Pulse4Power = 692.13*((0.351/Wavelength)^2)*SizeMod^2;
    RiseTime = 2e-10;
    LaserEndTime = [11].*10^-9;
    
    %Here I am changing the vapour density by giving it an array (while not
    %changing any other variable)
    RadiusVapour = [0.125];
    RadiusIce = [0.137];
    RadiusCH=0.1425;
    DensityVapour=[0.0012:0.0001:0.0020];
    
    ChangeTimeAt = [0 5 8 11 15].*10^-9;
    TimeSpacing = [1E-10 1E-11 1E-10 1E-11];
    
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
    SecondLaser, SecondLaserOnTime, LaserEndTime, SecondLaserWavelength, SecondLaserPower, Wavelength);
end


%% Example 8 - Changing LaserEndTime.
if Selection==108
    
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
    
    %This time we vary the LaserEndTime. Both the main laser and the second
    %laser will begin to switch off at this time
    Pulse1Power = 2*((0.351/Wavelength)^2)*SizeMod^2;
    Pulse2Power = 14*((0.351/Wavelength)^2)*SizeMod^2;;
    Pulse3Power = 98*((0.351/Wavelength)^2)*SizeMod^2;
    Pulse4Power = 692.13*((0.351/Wavelength)^2)*SizeMod^2;
    RiseTime = 2e-10;
    LaserEndTime = [11:0.2:15].*10^-9;
    
    RadiusVapour = [0.125];
    RadiusIce = [0.137];
    RadiusCH=0.1425;
    DensityVapour=[0.0012];
    
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


%% Simple functions

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
Analysis_Git
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


%% The main functions
function [NumLayers, MaxMassDiff, Method, Regions, ...
    Time, PulseOnTime, PulsePower] = GeneratorFunction(...
    RadiusVapour, RadiusIce, RadiusCH, DensityVapour, FoamOnly, TNBurn,...
    Pulse1Time, Pulse2Time, Pulse3Time, Pulse4Time, Pulse1Power, Pulse2Power, Pulse3Power, Pulse4Power, ...
    RelativeTimingsMode, RiseTime, LaserEndTime, ChangeTimeAt, TimeSpacing,...
    ElectronHeatingOn, HeatingTimings, Heating, HeatingZones, HeatingEnergy, ...
     SecondLaser, SecondLaserOnTime, SecondLaserOffTime, SecondLaserWavelength, SecondLaserPower, Wavelength)


%Mode 1 - Thickness Scan
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
    
    
%Mode 2 - Laser Scan with 3 pulses, or 4 pulses changing 2nd and 3rd
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
    
    
    %Mode 3 - Laser Scan changing 3rd and 4th
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
    
%Mode 4 - Laser Scan with 4 pulses, changing 3rd and 4th pulses
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
    
    
    
%Mode=5 - Changing Electron Heating
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
    
%Mode=6 - ChangingLaserEndTime
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
    
    
    %Mode=7 - ChangingVapourPressure
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
    
    
        %Mode=8 - Should work if only second laser on time is being changed
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
