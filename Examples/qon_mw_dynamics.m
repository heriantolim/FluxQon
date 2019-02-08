%% FluxQon Example #4: Microwave Dynamics of Flux Qon
% 
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 04/11/2018
% Last modified: 06/11/2018

% Add the required packages using MatVerCon.
% addpackage('MatCommon','MatGraphics','PhysConst','FluxQon');

% Clear workspace variables.
clear;

NUM_PHOTONS=2;
QUBIT_DETUNING=.5;
ION_CONCENTRATION=.005;

% Create a qubit object.
qb=Circle3JJ(...
	'Frequency',2*pi*1.95e9*(1+QUBIT_DETUNING),...
	'TunnelingFrequency',2*pi*1014.7424e9,...
	'Radius',5e-6,...% 5 um
	'DecayRate',2*pi*1e5,...% 100 MHz
	'DephasingRate',2*pi*2e5,...% 200 MHz
	'Temperature',2e-2);% 20 mK

% Create a YSOEr object.
ion=CylinderYSOEr(...
	'FockDimension',NUM_PHOTONS+2,...
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
	'FockDimension',NUM_PHOTONS+2,...
	'DecayRate',2*pi*1e4,...% 10 kHz
	'Temperature',2e-2);% 20 mK
mw.MagneticAmplitude=ion.CouplingStrength*sqrt(ion.NumIons) ...
	*Constant.FluxQuantum/(pi*qb.Area*qb.TunnelingFrequency);

% Get the Hilbert dimensions of each object.
d=Hilbert.dimension(qb,ion,mw);
M=numel(d);

% Construct the initial state.
psi0=cell(2,M);
for i=1:M
	psi0{1,i}=i;% the state index.
	psi0{2,i}=eye(d(i),1);% the ground state.
end
% Excite the microwave state.
psi0{2,3}([1,NUM_PHOTONS+1])=psi0{2,3}([NUM_PHOTONS+1,1]);
psi0_mw=psi0{2,3};% the initial state of the microwave.
psi0=State.kron(d,psi0{:});% the initial state of the system.

% The observables.
O=cell(1,4);
O{1}=Operator.kron(d,1,qb.Number);
O{2}=Operator.kron(d,2,ion.Number(2));
O{3}=Operator.kron(d,3,mw.Number);
O{4}=@(x)fidelity(psi0_mw,subdsmat(x,d,2));

% Solve the Lindblad Master Equation. This might take a few minutes.
[time,obs]=Solve.LME(...
	5e-7,...% time span in seconds.
	psi0*psi0',...% the initial density matrix of the system.
	Construct.Lindblad(qb,ion,mw),...% the Lindblad superoperator.
	'Observable',O,'DispProgress',true);
time=time*1e9;% convert second to nanosecond.

%% Plotting
% Settings.
Groot.usedefault();
Groot.usedefault('latex',8,.6);
RESOLUTION=300;

L=Layout.tiledaxes(1,2,[5,3],[.1,.1,.1,.1],[-.3,.6,.6,.4],[.3,.3,.3,.3]);

X_LIM=[0,500];
Y_LIM={[-.1,1.3].*NUM_PHOTONS,[-.1,1.3]};
X_TICK=0:100:400;
LEGEND_STRING={...
	'$\langle\hat{\sigma}_+\hat{\sigma}_-\rangle$',...
	'$\langle\hat{b}_1^\dagger\hat{b}_1\rangle$',...
	'$\langle\hat{a}_\mathrm{mw}^\dagger\hat{a}_\mathrm{mw}\rangle$'};
PLOT_LINE_COLOR={...
	[40,80,200]/255,...% blue
	[20,200,40]/255,...% green
	[125,50,150]/255,...% purple
	[0,120,190]/255};% blue

% Figure.
fig=docfigure(L.Paper.Position(3:4));

% Axes.
bgaxes('Position',L.Container.Position);
xlabel('Time (ns)');

ax=cell(1,2);
for j=1:2
	ax{j}=axes('Position',L.Axes.Position{j},...
		'XLim',X_LIM,'YLim',Y_LIM{j},'XTick',X_TICK);
	fixticklength(.2);
end
ax{2}.YAxisLocation='right';
ylabel(ax{1},'Excitation number');
ylabel(ax{2},'Fidelity');

% Plots.
for k=1:3
	plot(ax{1},time,obs{k},'Color',PLOT_LINE_COLOR{k});
end
plot(ax{2},time,obs{4},'Color',PLOT_LINE_COLOR{4});

% Legend.
h=legend(ax{1},LEGEND_STRING,'Location','NorthEast');

% Figure labels.
Label.subfigure(ax{2},'Transfer fidelity','Position',.4);

% Saving.
print(fig,'qon_mw_dynamics.png','-dpng',sprintf('-r%d',RESOLUTION));
close(fig);
