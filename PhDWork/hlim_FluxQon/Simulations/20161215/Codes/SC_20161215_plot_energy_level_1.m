function SC_20161215_plot_energy_level_1(varargin)
%% Settings
Groot.usedefault();
Groot.usedefault('latex',12,.72);

PAPER_POSITION=[0,0,15.4,6.2];
AXES_POSITION={[2.1,1.3,8,4.8],[10.5,1.3,4.8,4.8]};
CONTAINER_POSITION=[.65,.7,14.8,5.9];
[M,N]=size(AXES_POSITION);

X_LIM={[-2.5,2.5],[-1.5,1.5]};
Y_LIM={[-1.5,1.5],[-1.5,1.5]};
TICK_LENGTH=.2;
X_TICK=-2:2;
Y_TICK=-1:.5:1;
X_TICK_LABEL={'-4$\,\mathrm{g}\sqrt{n}$','-2$\,\mathrm{g}\sqrt{n}$','0',...
	'2$\,\mathrm{g}\sqrt{n}$','4$\,\mathrm{g}\sqrt{n}$'};
Y_TICK_LABEL={'-2$\,\hbar\mathrm{g}\sqrt{n}$','-$\hbar\mathrm{g}\sqrt{n}$',...
	'0','$\hbar\mathrm{g}\sqrt{n}$','2$\,\hbar\mathrm{g}\sqrt{n}$'};
X_LABEL_STRING='Detuning, $\Delta$';
Y_LABEL_STRING='Energy';

FIGURE_LABEL_LOCATION={'North','NorthWest'};
FIGURE_LABEL_STRING={'(a)','(b)'};

LEGEND_LOCATION='East';
LEGEND_STRING={'$E_{n,+}-n\hbar\omega_c$','$E_{n,-}-n\hbar\omega_c$'};

PLOT_LINE_COLOR={[220,80,20]/255,[0,120,190]/255,[.3,.3,.3]};
PLOT_LINE_STYLE={'-','-','--'};

ANNOTATION_COLOR=[30,160,20]/255;
ANNOTATION_LINE_STYLE='--';

%% Input Parsing
P=inputParser;
P.addParameter('Rewrite','yes',@(x)any(strcmpi(x,{'yes','no','ask'})));
P.addParameter('FileFormat',{'.fig','.pdf'},@FileDir.isfileext);
P.addParameter('Resolution',[300,600],@isintegervector);
P.parse(varargin{:});
P=P.Results;

%% Data Processing
X=linspace(X_LIM{1}(1),X_LIM{1}(2),ceil(300/2.54*AXES_POSITION{1}(3)));
Y=cell(N,2);
Y{1,1}=+.5*sqrt(1+X.^2);
Y{1,2}=-.5*sqrt(1+X.^2);
Y{2,1}=Y{1,1}+.5*X;
Y{2,2}=Y{1,2}+.5*X;

%% Drawing
fig=docfigure(PAPER_POSITION(3:4));

bgaxes('Position',CONTAINER_POSITION);
xlabel(X_LABEL_STRING);
ylabel(Y_LABEL_STRING);

ax=cell(M,N);
for i=1:M
	for j=1:N
		ax{i,j}=axes('Position',AXES_POSITION{i,j},...
			'XLim',X_LIM{i,j},'YLim',Y_LIM{i,j},...
			'XTick',X_TICK,'YTick',Y_TICK,...
			'XTickLabel',X_TICK_LABEL,'YTickLabel',Y_TICK_LABEL);
		fixticklength(TICK_LENGTH);
		if i~=M
			ax{i,j}.XTickLabel='';
		end
		if j~=1
			ax{i,j}.YTickLabel='';
		end
	end
end

plot(ax{1},X_LIM{1},X_LIM{1}.'*[-.5,.5],...
	'Color',PLOT_LINE_COLOR{3},'LineStyle',PLOT_LINE_STYLE{3});
plot(ax{2},X_LIM{2},[X_LIM{2};0,0].',...
	'Color',PLOT_LINE_COLOR{3},'LineStyle',PLOT_LINE_STYLE{3});

pl=cell(N,2);
for j=1:N
	for k=1:2
		pl{j,k}=plot(ax{j},X,Y{j,k},...
			'Color',PLOT_LINE_COLOR{k},'LineStyle',PLOT_LINE_STYLE{k});
	end
end

legend([pl{1,:}],LEGEND_STRING,'Location',LEGEND_LOCATION);

%% Annotating
[px,py]=data2pos(ax{1},[X_LIM{1}(1),0],[-.5,.5]);
for k=1:2
	annotation('line','Position',[px(1),py(k),px(2)-px(1),0],...
		'Color',ANNOTATION_COLOR,'LineStyle',ANNOTATION_LINE_STYLE);
end

for i=1:M
	for j=1:N
		Label.subfigure(ax{i,j},FIGURE_LABEL_STRING{i,j},...
			'Location',FIGURE_LABEL_LOCATION{i,j},'Position',.4);
	end
end

%% Saving
savefigure(fig,'F','energy','level',1,...
	P.FileFormat,'Resolution',P.Resolution,'Rewrite',P.Rewrite);
close(fig);

end
