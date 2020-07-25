function SC_20161215_plot_Rabi_oscillation_1(varargin)
%% Settings
Groot.usedefault();
Groot.usedefault('latex',12,.72);

PAPER_POSITION=[0,0,14.9,7];
AXES_POSITION={[0,0,7,7],[9.2,1.3,5.5,5.5]};
[M,N]=size(AXES_POSITION);

X_LIM=[0,9];
Y_LIM=[0,1];
TICK_LENGTH=.2;
X_TICK=0:2:8;
Y_TICK=.25:.25:1;
X_TICK_LABEL={'0','','8$\,\mathrm{g}\sqrt{n}$','','16$\,\mathrm{g}\sqrt{n}$'};
Y_TICK_LABEL={'$\frac{1}{4}\sin\theta$','$\frac{1}{2}\sin\theta$',...
	'$\frac{3}{4}\sin\theta$','$\sin\theta$'};
X_LABEL_STRING='Detuning, $|\Delta|$';
Y_LABEL_STRING='Amplitude';

FIGURE_LABEL_LOCATION={'NorthWest','NorthEast'};
FIGURE_LABEL_STRING={'(a)','(b)'};

PLOT_LINE_COLOR=[0,120,190]/255;

%% Input Parsing
P=inputParser;
P.addParameter('Rewrite','yes',@(x)any(strcmpi(x,{'yes','no','ask'})));
P.addParameter('FileFormat',{'.fig','.pdf'},@FileDir.isfileext);
P.addParameter('Resolution',[300,600],@isintegervector);
P.parse(varargin{:});
P=P.Results;

%% Data Processing
X=linspace(X_LIM(1),X_LIM(2),ceil(300/2.54*AXES_POSITION{2}(3)));
O=sqrt(1+X.^2);
Y=sin(2*atan(sqrt((O+X)./(O-X))));

%% Drawing
fig=docfigure(PAPER_POSITION(3:4));

ax=cell(M,N);
ax{1}=axes('Position',AXES_POSITION{1},...
	'XLim',[0,1],'YLim',[0,1],...
	'XTick',[],'YTick',[],...
	'XColor','w','YColor','w');

ax{2}=axes('Position',AXES_POSITION{2},...
	'XLim',X_LIM,'YLim',Y_LIM,...
	'XTick',X_TICK,'YTick',Y_TICK,...
	'XTickLabel',X_TICK_LABEL,'YTickLabel',Y_TICK_LABEL);
fixticklength(TICK_LENGTH);
xlabel(X_LABEL_STRING);
ylabel(Y_LABEL_STRING);

axes(ax{1});
imagesc([0,1],[1,0],imread(fullpath('F',...
	'Rabi','Bloch','sphere','1','600dpi','.png')));

plot(ax{2},X,Y,'Color',PLOT_LINE_COLOR);

for i=1:M
	for j=1:N
		Label.subfigure(ax{i,j},FIGURE_LABEL_STRING{i,j},...
			'Location',FIGURE_LABEL_LOCATION{i,j},'Position',.4);
	end
end

%% Saving
savefigure(fig,'F','Rabi','oscillation',1,...
	P.FileFormat,'Resolution',P.Resolution,'Rewrite',P.Rewrite);
close(fig);

end