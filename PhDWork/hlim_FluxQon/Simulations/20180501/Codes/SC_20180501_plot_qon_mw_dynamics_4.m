function SC_20180501_plot_qon_mw_dynamics_4(n,varargin)
%% Settings
Groot.usedefault();
Groot.usedefault('latex',12,.72);

M=2; N=2;
L=Layout.tiledaxes(M,N,[6.4,4],[.1,.1,.1,.1],[.4,.7,.7,.7],[.4,.4,.4,.4]);
fprintf('Paper width  = %g\n',L.Paper.Position(3));
fprintf('Paper height = %g\n',L.Paper.Position(4));

X_LIM=[0,500];
Y_LIM={[-.1*n,1.3*n],[-.1,1.3]};
TICK_LENGTH=.2;
X_TICK=0:100:400;
X_LABEL_STRING='Time (ns)';
Y_LABEL_STRING={'Excitation number','Fidelity'};
FIGURE_LABEL_STRING={...
	'(a) $\Delta=\frac{1}{2}\,\omega$','(b)';
	'(c) $\Delta=2\,\omega$','(d)'};

LEGEND_STRING={{...
	'$\langle\hat{\sigma}_+\hat{\sigma}_-\rangle\quad$',...
	'$\langle\hat{b}_1^\dagger\hat{b}_1\rangle\quad$',...
	'$\langle\hat{a}_\mathrm{mw}^\dagger\hat{a}_\mathrm{mw}\rangle$'},{...
	'Transfer fidelity'}};

PLOT_LINE_COLOR={{...
	[40,80,200]/255,...% blue
	[20,200,40]/255,...% green
	[125,50,150]/255},{...% purple
	[0,120,190]/255,...% blue
	[220,80,20]/255}};% red

%% Input Parsing
P=inputParser;
P.addRequired('NumPhotons',@isintegerscalar);
P.addOptional('CouplingRatio',1,@isrealscalar);
P.addOptional('Concentration',.005,@isrealscalar);% atomic percent
P.addParameter('Rewrite','yes',@(x)any(strcmpi(x,{'yes','no','ask'})));
P.addParameter('FileFormat',{'.fig','.pdf'},@FileDir.isfileext);
P.addParameter('Resolution',[300,600],@isintegervector);
P.parse(n,varargin{:});
P=P.Results;

%% Data Processing
D=[.5,2];
X=cell(M,1);
Y=cell(M,N);
for i=1:M
	S=loaddata('R','qon','mw','dynamics',n,2,D(i),...
		P.CouplingRatio,P.Concentration);
	X{i}=S.Time*1e9;
	Y{i,1}=S.Number;
	Y{i,2}={S.Fidelity};
end

%% Drawing
fig=docfigure(L.Paper.Position(3:4));

bgaxes('Position',L.Container.Position);
xlabel(X_LABEL_STRING);
ylabel(Y_LABEL_STRING{1});
bgaxes('Position',L.Container.Position,'YAxisLocation','right');
ylabel(Y_LABEL_STRING{2});

ax=cell(M,N);
for i=1:M
	for j=1:N
		ax{i,j}=axes('Position',L.Axes.Position{i,j},...
			'XLim',X_LIM,'YLim',Y_LIM{j},'XTick',X_TICK);
		fixticklength(TICK_LENGTH);
		if i~=M
			ax{i,j}.XTickLabel='';
		end
		if j~=1
			ax{i,j}.YAxisLocation='right';
		end
		for k=1:numel(Y{i,j})
			plot(X{i},Y{i,j}{k},'Color',PLOT_LINE_COLOR{j}{k});
		end
	end
end

h=legend(ax{1,1},LEGEND_STRING{1},'NumColumns',3);
p=h.Position;
p(1)=ax{1,1}.Position(1)+ax{1,1}.Position(3)/2-p(3)/2;
p(2)=sum(ax{1,1}.Position([2,4]))+.2;
h.Position=p;

h=legend(ax{1,2},LEGEND_STRING{2});
p=h.Position;
p(1)=ax{1,2}.Position(1)+ax{1,2}.Position(3)/2-p(3)/2;
p(2)=sum(ax{1,2}.Position([2,4]))+.2;
h.Position=p;

%% Annotating
for i=1:M
	for j=1:N
		Label.subfigure(ax{i,j},FIGURE_LABEL_STRING{i,j},'Position',.4);
	end
end

%% Saving
savefigure(fig,'F','qon','mw','dynamics',4,n,...
	P.FileFormat,'Resolution',P.Resolution,'Rewrite',P.Rewrite);
close(fig);

end