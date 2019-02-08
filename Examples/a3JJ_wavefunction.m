%% FluxQon Example #2: 3JJ Qubit Wavefunction
% 
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 04/11/2018
% Last modified: 04/11/2018

% Add the required packages using MatVerCon.
% addpackage('MatCommon','MatGraphics','PhysConst','FluxQon');

% Clear workspace variables.
clear;

% Create the 3JJ flux qubit object.
JosephsonEnergy=Constant.Planck*785.2552e9;
CoulombEnergy=JosephsonEnergy/10;
qb=a3JJ(CoulombEnergy,JosephsonEnergy,.7);

% Solve the time-independent Schrodinger equation on 100x100 coordinates.
% This calculation might take a few minutes.
[~,W]=qb.solveTISE(100);
X=linspace(-pi/2,pi/2,size(W,1));

% Print the qubit frequencies to console.
fprintf('Excitation Frequency = %.6f GHz\n',qb.Frequency/2/pi/1e9);
fprintf('Tunneling Frequency  = %.6f GHz\n',qb.TunnelingFrequency/2/pi/1e9);

%% Plotting
% Settings.
Groot.usedefault();
Groot.usedefault('latex',8,.6);
RESOLUTION=300;

L=Layout.tiledaxes(1,2,[5,5],[.1,.4,1.2,1.6],[0,0,-.3,-.5]);

X_LIM=[-pi/2,pi/2];
X_TICK=-pi/2:pi/4:pi/2;
X_TICK_LABEL={'-$\frac{\pi}{2}$','-$\frac{\pi}{4}$','0',...
	'$\frac{\pi}{4}$','$\frac{\pi}{2}$'};
Z_LIM=[min(min(min(W(:,:,1))),min(min(W(:,:,2)))),...
	max(max(max(W(:,:,1))),max(max(W(:,:,2))))]*[1.25,-.25;-.25,1.25];
FIGURE_LABEL_STRING={'Ground state','Excited state'};

% Figure.
fig=docfigure(L.Paper.Position(3:4));

% Axes.
ax=cell(1,2);
for j=1:2
	ax{j}=axes('Position',L.Axes.Position{j},'ZLim',Z_LIM,'ZTick',-2:2,...
		'XLim',X_LIM,'XTick',X_TICK,'XTickLabel',X_TICK_LABEL,...
		'YLim',X_LIM,'YTick',X_TICK,'YTickLabel',X_TICK_LABEL);
	fixticklength(.2);
	grid(ax{j},'on');
	view(ax{j},[1,5,1]);
	xlabel('$\phi_+$');
	ylabel('$\phi_-$');
	zlabel('Wavefunction, $\psi$');
end

% Plots.
mesh(ax{2},X,X,W(:,:,2).');
ax{1}.CLim=ax{2}.CLim;
mesh(ax{1},X,X,W(:,:,1).');

% Figure labels.
for j=1:2
	Label.subfigure(ax{j},FIGURE_LABEL_STRING{j},'Position',.4);
end

% Saving.
print(fig,'a3JJ_wavefunction.png','-dpng',sprintf('-r%d',RESOLUTION));
close(fig);
