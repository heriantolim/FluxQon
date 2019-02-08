%% FluxQon Example #1: YSO:Er Energy vs Angle
% 
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 04/11/2018
% Last modified: 04/11/2018

% Add the required packages using MatVerCon.
% addpackage('MatCommon','MatGraphics','PhysConst','FluxQon');

% Clear workspace variables.
clear;

% Create a YSOEr object with magnetic class orientation set to 1
% and bias magnetic field set to 80 mT.
ion=YSOEr();
ion.Class=1;
ion.MagneticField=.08;

% Calculate the energy of the YSOEr for isotope 1 and 2, site 1 and 2, and
% for various Euler rotation angles.
N=181;
angle=linspace(0,pi,N);
energy=cell(2,2);
for i=1:2
	for j=1:2
		ion.Isotope=i;
		ion.Site=j;
		energy{i,j}=zeros(ion.HilbertDimension,N);
		for n=1:N
			ion.Rotation=[angle(n),pi/2,0];
			energy{i,j}(:,n)=diag(ion.Hamiltonian)/Constant.Planck/1e9;% in GHz
		end
	end
end
angle=angle*180/pi;

%% Plotting
% Settings.
Groot.usedefault();
Groot.usedefault('latex',8,.6);
RESOLUTION=300;

L=Layout.tiledaxes(1,2,[5,5],[.1,.1,.3,.3],[-.2,-.3,.5,.5],[.3,.3,.1,.1]);

FIGURE_LABEL_STRING={'Site 1','Site 2'};
PLOT_LINE_COLOR={...
	{[228, 26, 28]/255,[ 55,126,184]/255,[ 77,175, 74]/255,[152, 78,163]/255,...
	 [255,127,  0]/255,[166, 86, 40]/255,[247,129,191]/255,[153,153,153]/255},...
	{'k'}};
PLOT_LINE_WIDTH=[.5,.7];

% Figure.
fig=docfigure(L.Paper.Position(3:4));

% Axes.
bgaxes('Position',L.Container.Position);
xlabel('$\gamma$ angle (degree)');
ylabel('Energy (GHz)');
ax=cell(1,2);
for j=1:2
	ax{j}=axes('Position',L.Axes.Position{j},'XLim',[angle(1),angle(end)]);
	fixticklength(.2);
end
ax{2}.YTickLabel='';

% Plots.
yLim=zeros(2,2);
for j=1:2
	for i=1:2
		for l=1:size(energy{i,j},1)
			plot(ax{j},angle,energy{i,j}(l,:),...
				'Color',PLOT_LINE_COLOR{i}{ceil(l/2)},...
				'LineWidth',PLOT_LINE_WIDTH(i));
		end
	end
	yLim(j,:)=ax{j}.YLim;
end

% Reconfigure the axes.
yLim=[min(yLim(:,1)),max(yLim(:,2))];
for j=1:2
	ax{j}.YLim=yLim;
end

% Figure labels.
for j=1:2
	Label.subfigure(ax{j},FIGURE_LABEL_STRING{j},'Position',.4);
end

% Saving.
print(fig,'YSOEr_energy_vs_angle.png','-dpng',sprintf('-r%d',RESOLUTION));
close(fig);
