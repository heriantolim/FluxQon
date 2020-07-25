function SC_20180502_plot_qon_dynamics_1(n,varargin)
%% Settings
Groot.usedefault();
Groot.usedefault('latex',12,.72);

M=1; N=2;
L=Layout.tiledaxes(M,N,[6.4,4],[.1,.1,.1,.1],[.4,-.4,.7,.7],[.4,.4,.4,.4]);
fprintf('Paper width  = %g\n',L.Paper.Position(3));
fprintf('Paper height = %g\n',L.Paper.Position(4));

X_LIM=[0,500];
Y_LIM=[-.1*n,1.3*n];
TICK_LENGTH=.2;
X_TICK=0:100:400;
X_LABEL_STRING='Time (ns)';
Y_LABEL_STRING='Excitation number';
FIGURE_LABEL_STRING={...
	'(a) Initial exc. in Er$^\mathrm{3+}$ Site 1',...
	'(b) Initial exc. in Er$^\mathrm{3+}$ Site 2'};

LEGEND_STRING={...
	'$\langle\hat{\sigma}_+\hat{\sigma}_-\rangle\quad$',...
	'$\langle\hat{b}_1^\dagger\hat{b}_1\rangle\quad$',....
	'$\langle\hat{c}_1^\dagger\hat{c}_1\rangle\quad$'};

PLOT_LINE_COLOR={...
	[40,80,200]/255,...% blue
	[20,200,40]/255,...% green
	[255,170,0]/255};% orange

%% Input Parsing
P=inputParser;
P.addRequired('NumPhotons',@isintegerscalar);
P.addOptional('Concentration',.005,@isrealscalar);% atomic percent
P.addParameter('Rewrite','yes',@(x)any(strcmpi(x,{'yes','no','ask'})));
P.addParameter('FileFormat',{'.fig','.pdf'},@FileDir.isfileext);
P.addParameter('Resolution',[300,600],@isintegervector);
P.parse(n,varargin{:});
P=P.Results;

%% Data Processing
F(1,1)=FileDir('R','qon','dynamics',n,2,1,P.Concentration);
F(1,2)=FileDir('R','qon','dynamics',n,2,2,P.Concentration);

X=cell(M,N);
Y=cell(M,N);
for i=1:M
	for j=1:N
		S=loaddata(F(i,j));
		X{i,j}=S.Time*1e9;
		Y{i,j}=S.Number;
	end
end

%% Drawing
fig=docfigure(L.Paper.Position(3:4));

bgaxes('Position',L.Container.Position);
xlabel(X_LABEL_STRING);
ylabel(Y_LABEL_STRING);

ax=cell(M,N);
for i=1:M
	for j=1:N
		ax{i,j}=axes('Position',L.Axes.Position{i,j},...
			'XLim',X_LIM,'YLim',Y_LIM,'XTick',X_TICK);
		fixticklength(TICK_LENGTH);
		if i~=M
			ax{i,j}.XTickLabel='';
		end
		if j~=1
			ax{i,j}.YTickLabel='';
		end
		for k=1:numel(Y{i,j})
			plot(X{i,j},Y{i,j}{k},'Color',PLOT_LINE_COLOR{k});
		end
	end
end
ax{2}.Children(2).LineStyle='--';

h=legend(ax{1},LEGEND_STRING,'NumColumns',3);
p=h.Position;
p(1)=mean([ax{1}.Position(1),sum(ax{2}.Position([1,3]))])-p(3)/2;
p(2)=max(sum(ax{1}.Position([2,4])),sum(ax{2}.Position([2,4])))+.2;
h.Position=p;

%% Annotating
for i=1:M
	for j=1:N
		Label.subfigure(ax{i,j},FIGURE_LABEL_STRING{i,j},'Position',.4);
	end
end

%% Saving
savefigure(fig,'F','qon','dynamics',1,n,...
	P.FileFormat,'Resolution',P.Resolution,'Rewrite',P.Rewrite);
close(fig);

end