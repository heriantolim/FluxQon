%% FluxQon Example #3: Unitary Down-conversion
% 
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 06/11/2018
% Last modified: 06/11/2018

% Add the required packages using MatVerCon.
% addpackage('MatCommon','MatGraphics','PhysConst','FluxQon');

% Clear workspace variables.
clear;

NUM_PHOTONS=1;
NUM_IONS=1e4;
QUBIT_DETUNING=0;
LINE_STRENGTH3=1e3;
COUPLING_RATIO=sqrt(NUM_IONS/NUM_PHOTONS);

% Create a qubit object.
qb=Circle3JJ(...
	'Frequency',2*pi*1.95e9*(1+QUBIT_DETUNING),...
	'TunnelingFrequency',2*pi*1014.7424e9,...
	'Radius',5e-6);% 5 um

% Create a YSOEr object.
ion=CylinderYSOEr(...
	'NumIons',NUM_IONS,...
	'FockDimension',NUM_PHOTONS+2,...
	'Isotope',2,...
	'Multiplet',1:2,...
	'Site',1,...
	'Class',1,...
	'Rotation',[152,270,0]*pi/180,...
	'MagneticField',0.1301185,...% 130 mT
	'Radius',4e-6,...% 4 um
	'Height',8e-6,...% 8 um
	'CouplingStrength',2*pi*[5645.1756,COUPLING_RATIO*LINE_STRENGTH3],...
	'LineStrength',2*pi*LINE_STRENGTH3*[0,10,1,1,10,0]);

% Create a microwave resonator object.
mw=Microwave(...
	'Frequency',2*pi*1.95e9,...
	'FockDimension',NUM_PHOTONS+2);
mw.MagneticAmplitude=2*pi*COUPLING_RATIO*LINE_STRENGTH3 ...
	*Constant.FluxQuantum/(pi*qb.Area*qb.TunnelingFrequency);

% Create an optical resonator object.
op=Optical(...
	'Frequency',diff(ion.Frequency([1,4])),...
	'FockDimension',NUM_PHOTONS+2);

% Get the Hilbert dimensions of each object.
d=Hilbert.dimension(qb,ion,mw,op);
M=numel(d);

% Construct the initial state.
psi0=cell(2,M);
for i=1:M
	psi0{1,i}=i;% the state index.
	psi0{2,i}=eye(d(i),1);% the ground state.
end
% Excite the optical state.
psi0{2,4}([1,NUM_PHOTONS+1])=psi0{2,4}([NUM_PHOTONS+1,1]);
psi0_op=psi0{2,4};% the initial state of the optical.
psi0=State.kron(d,psi0{:});% the initial state of the system.

% The observables.
O=cell(1,3);
O{1}=Operator.kron(d,3,mw.Number);
O{2}=Operator.kron(d,4,op.Number);
O{3}=@(x)fidelity(psi0_op,subdsmat(x,d,3));

% Solve the Schrodinger equation. This might take a few minutes.
[time,obs]=Solve.UTE(...
	0:1e-7:1e-4,...% time points in seconds.
	psi0,...% the initial state of the system.
	Construct.Hamiltonian(qb,ion,mw,op),...% the Hamiltonian of the system.
	'Observable',O,'DispProgress',true);
time=time*1e6;% convert second to microsecond.

%% Plotting
% Settings.
Groot.usedefault();
Groot.usedefault('latex',8,.6);
RESOLUTION=300;

L=Layout.tiledaxes(1,2,[5,3],[.1,.1,.1,.1],[-.3,.6,.6,.6],[.3,.3,.3,.3]);

X_LIM=[0,100];
Y_LIM={[-.1,1.3].*NUM_PHOTONS,[-.1,1.3]};
X_TICK=0:20:80;
LEGEND_STRING={...
	'$\langle\hat{a}_\mathrm{mw}^\dagger\hat{a}_\mathrm{mw}\rangle$',...
	'$\langle\hat{a}_\mathrm{op}^\dagger\hat{a}_\mathrm{op}\rangle$'};
PLOT_LINE_COLOR={...
	[125,50,150]/255,...% purple
	[255,130,40]/255,...% orange
	[0,120,190]/255};% blue

% Figure.
fig=docfigure(L.Paper.Position(3:4));

% Axes.
bgaxes('Position',L.Container.Position);
xlabel('Time ($\mu$s)');

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
for k=1:2
	plot(ax{1},time,obs{k},'Color',PLOT_LINE_COLOR{k});
end
plot(ax{2},time,obs{3},'Color',PLOT_LINE_COLOR{3});

% Legend.
h=legend(ax{1},LEGEND_STRING,'Location','NorthEast','NumColumns',2);

% Figure labels.
Label.subfigure(ax{2},'Down-conversion fidelity','Position',.4);

% Saving.
print(fig,'unitary_downconversion.png','-dpng',sprintf('-r%d',RESOLUTION));
close(fig);
