%% FluxQon Example #5: Transmission Spectrum
% 
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 06/11/2018
% Last modified: 06/11/2018

% Add the required packages using MatVerCon.
% addpackage('MatCommon','MatGraphics','PhysConst','FluxQon');

% Clear workspace variables.
clear;

ION_CONCENTRATION=.005;

% Create a qubit object.
qb=Circle3JJ(...
	'Frequency',2*pi*1.95e9,...
	'TunnelingFrequency',2*pi*1014.7424e9,...
	'Radius',5e-6,...% 5 um
	'DecayRate',2*pi*1e5,...% 100 MHz
	'DephasingRate',2*pi*2e5,...% 200 MHz
	'Temperature',2e-2);% 20 mK

% Create a YSOEr object.
ion=CylinderYSOEr(...
	'FockDimension',2,...
	'Isotope',2,...
	'Multiplet',1,...
	'Site',1,...
	'Class',1,...
	'Rotation',[152,270,0]*pi/180,...
	'MagneticField',0.0835351,...% 83.5 mT
	'Radius',4e-6,...% 4 um
	'Height',8e-6,...% 8 um
	'CouplingStrength',2*pi*5645.1755,...
	'LineStrength',2*pi*10,...% 10 Hz
	'DecayRate',2*pi*1e6,...% 1 MHz
	'Temperature',2e-2);% 20 mK
% NumIons = 1/2 * Concentration * Volume * NumberDensity
% half because there are two types of Er ions in the crystal.
ion.NumIons=round(ION_CONCENTRATION*ion.Volume/200 ...
	*Constant.AvogadroNumber*ion.Density/ion.MolarMass);

% Create a microwave resonator object.
mw=Microwave(...
	'Frequency',2*pi*1.95e9,...
	'FockDimension',2,...
	'DecayRate',2*pi*1e4,...% 10 kHz
	'Temperature',2e-2);% 20 mK
mw.MagneticAmplitude=ion.CouplingStrength*sqrt(ion.NumIons) ...
	*Constant.FluxQuantum/(pi*qb.Area*qb.TunnelingFrequency);

% Calculating the transmission/radiation spectrum.
% This might take a few minutes.
S=RadiationSpectrum(...
	2*1.95e9,...% maximum probe frequency in Hz.
	1e-4,...% probe time length in second.
	Construct.Lindblad(qb,ion,mw),...% the Lindblad superoperator.
	Hilbert.dimension(qb,ion,mw),...% the Hilbert dimensions.
	3,...% set the readout operator to be the annihilation operator of mw object.
	'DispProgress',true);
S(1,:)=S(1,:)/1.95e9;% convert the units of the frequency.

%% Plotting
% Settings.
Groot.usedefault();
Groot.usedefault('latex',8,.6);
RESOLUTION=300;
AXES_SIZE=[12,4];
X_LIM=[.9,1.1];
Y_LIM=[0,2e-5];
X_TICK=.9:.05:1.1;
X_TICK_LABEL=arrayfun(@(x)sprintf('%.2f$\\,\\omega_0$',x),...
	X_TICK,'UniformOutput',false);
X_TICK_LABEL{3}='$\omega_0$';

% Figure.
fig=docfigure(AXES_SIZE);

% Axes.
pos=[0,0,AXES_SIZE];
ax=axes('Position',pos,'XLim',X_LIM,'YLim',Y_LIM,...
	'XTick',X_TICK,'XTickLabel',X_TICK_LABEL);
xlabel('Frequency');
ylabel('Intensity');
fixticklength(.2);

% Plots.
plot(S(1,:),S(2,:));

% Reconfigure the layout.
margin=ax.TightInset+.1;
ax.Position=pos+[margin(1),margin(2),0,0];
pos=pos+[0,0,margin(1)+margin(3),margin(2)+margin(4)];
set(fig,{'Position','PaperPosition','PaperSize'},{pos,pos,pos(3:4)});

% Saving.
print(fig,'transmission_spectrum.png','-dpng',sprintf('-r%d',RESOLUTION));
close(fig);
