function SC_20160625_plot_state_and_energy_1(varargin)
%% Settings
Groot.usedefault();
Groot.usedefault('latex',12,.72);

M=1; N=2; O=2;
L=Layout.tiledaxes(M,N,[5,5],[.5,1.7,1.5,1.5],[.1,-1.3,-.2,-.3]);
fprintf('Paper width  = %g\n',L.Paper.Position(3));
fprintf('Paper height = %g\n',L.Paper.Position(4));

X_LIM={[0,2*pi],[-2*pi,2*pi]};
Y_LIM={{[-1,2],[.5,2]}};
TICK_LENGTH=.2;
X_TICK={0:pi/2:2*pi,-2*pi:pi:2*pi};
X_TICK_LABEL={{'0','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','2$\pi$'},...
	{'-2$\pi$','-$\pi$','0','$\pi$','2$\pi$'}};
X_LABEL_STRING={'Phase shift, $\phi$',...
	'Bias phase, $(\phi_\mathrm{bias}-\pi)\times 10^3$'};
Y_LABEL_STRING={{'Wavefunction, $\psi(\phi)$','Potential, $U/E_J$'},...
	{'Energy, $(E-\bar{E})/E_J$'}};

LEGEND_LOCATION={'NorthEast','North'};
LEGEND_STRING={{'Potential','Ground','1\textsuperscript{st} excited'},...
	{'Ground','1\textsuperscript{st} excited'}};

FIGURE_LABEL_LOCATION={'NorthWest','West'};
FIGURE_LABEL_STRING={'(a)','(b)'};

PLOT_LINE_COLOR={'k',[0,120,190]/255,[220,80,20]/255};
PLOT_LINE_STYLE={'--','-','-'};

%% Input Parsing
P=inputParser;
P.addOptional('InductanceEnergy',2/pi,@isrealscalar);% time Josephson energy
P.addParameter('Rewrite','yes',@(x)any(strcmpi(x,{'yes','no','ask'})));
P.addParameter('FileFormat',{'.fig','.pdf'},@FileDir.isfileext);
P.addParameter('Resolution',[300,600],@isintegervector);
P.parse(varargin{:});
P=P.Results;

%% Data Processing
X=cell(M,N);
Y=X;

S=loaddata('R','state','at','pi');
X{1}=linspace(0,2*pi,size(S,1));
Y{1}=S(:,1:O).';
U=P.InductanceEnergy*(X{1}-pi).^2/2-cos(X{1});
U=diff(Y_LIM{1}{1})/diff(Y_LIM{1}{2})*(U-Y_LIM{1}{2}(1))+Y_LIM{1}{1}(1);

S=loaddata('R','energy','near','pi');
X{2}=(S(1,:)-pi)*1e3;
Y{2}=S(2:(O+1),:);
Y{2}=Y{2}-interp1(X{2},mean(Y{2},1),0);

%% Drawing
fig=docfigure(L.Paper.Position(3:4));

ax=cell(M,N);
for i=1:M
	for j=1:N
		ax{i,j}=axes('Position',L.Axes.Position{i,j},...
			'XLim',X_LIM{i,j},'XTick',X_TICK{i,j},'XTickLabel',X_TICK_LABEL{i,j});
		fixticklength(TICK_LENGTH);
		xlabel(X_LABEL_STRING{i,j});
		ylabel(Y_LABEL_STRING{i,j}{1});
	end
end
set(ax{1},'Box','off','YLim',Y_LIM{1}{1});
set(ax{2},'YLim',[min(Y{2}(:)),max(Y{2}(:))]);

plot(ax{1},X{1},U,'Color',PLOT_LINE_COLOR{1},'LineStyle',PLOT_LINE_STYLE{1});
for i=1:M
	for j=1:N
		for k=1:O
			plot(ax{i,j},X{i,j},Y{i,j}(k,:),...
				'Color',PLOT_LINE_COLOR{k+1},'LineStyle',PLOT_LINE_STYLE{k+1});
		end
		legend(ax{i,j},LEGEND_STRING{i,j},'Location',LEGEND_LOCATION{i,j});
	end
end

for i=1:M
	for j=1:N
		Label.subfigure(ax{i,j},FIGURE_LABEL_STRING{i,j},...
			'Location',FIGURE_LABEL_LOCATION{i,j},'Position',.4);
	end
end

axes('Position',L.Axes.Position{1},'Box','off','Color','none',...
	'XLim',X_LIM{1},'YLim',Y_LIM{1}{2},'XTick',X_TICK{1},'XTickLabel','',...
	'XAxisLocation','top','YAxisLocation','right');
fixticklength(TICK_LENGTH);
ylabel(Y_LABEL_STRING{1}{2});

%% Saving
savefigure(fig,'F','state','and','energy',1,...
	P.FileFormat,'Resolution',P.Resolution,'Rewrite',P.Rewrite);
close(fig);

end
