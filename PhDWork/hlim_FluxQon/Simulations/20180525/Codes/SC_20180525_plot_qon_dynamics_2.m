function SC_20180525_plot_qon_dynamics_2(varargin)
%% Settings
Groot.usedefault();
Groot.usedefault('latex',12,.72);

M=2; N=2;
L=Layout.tiledaxes(M,N,[6.4,4],[.1,.1,.1,.1],[.9,.7,.7,.7],[.4,.4,.4,.4]);
fprintf('Paper width  = %g\n',L.Paper.Position(3));
fprintf('Paper height = %g\n',L.Paper.Position(4));

X_LIM=[0,10];
Y_LIM=[-.1,1.3];
TICK_LENGTH=.2;
X_TICK=0:2:8;
X_LABEL_STRING='Time (ms)';
Y_LABEL_STRING={'Excitation number','Fidelity'};
FIGURE_LABEL_STRING={...
	'(a) $n=1$, $N=2$','(b)';
	'(c) $n=2$, $N=2$','(d)'};

LEGEND_STRING={{...
	'$\langle\hat{\sigma}_+\hat{\sigma}_-\rangle\quad$',...
	'$\langle\sum_i\hat{T}_{(i),3,3}\rangle\quad$',...
	'$\langle\hat{a}_\mathrm{mw}^\dagger\hat{a}_\mathrm{mw}\rangle$',...
	'$\langle\hat{a}_\mathrm{op}^\dagger\hat{a}_\mathrm{op}\rangle$'},{...
	'Down-conversion fidelity'}};

PLOT_LINE_COLOR={{...
	[40,80,200]/255,...% blue
	[20,225,130]/255,...% green
	[125,50,150]/255,...% purple
	[255,130,40]/255},{...% orange
	[0,120,190]/255}};% blue

%% Input Parsing
P=inputParser;
P.addOptional('CalcMode',0,@isintegerscalar);
P.addOptional('LineStrength3',1e3,@isrealscalar);
P.addParameter('Rewrite','yes',@(x)any(strcmpi(x,{'yes','no','ask'})));
P.addParameter('FileFormat',{'.fig','.pdf'},@FileDir.isfileext);
P.addParameter('Resolution',[300,600],@isintegervector);
P.parse(varargin{:});
P=P.Results;

%% Data Processing
nn=1:2;
NN=2;

X=cell(M,1);
Y=cell(M,N);
yLim=cell(M,N);
for i=1:2
	S=loaddata('R','qon','dynamics',P.CalcMode,...
		nn(i),NN,1,0,sqrt(NN/nn(i)),P.LineStrength3);
	X{i}=S.Time*1e3;
	Y{i,1}=[S.Number{1}; S.Number{5}+S.Number{9}; S.Number{10}; S.Number{11}];
	Y{i,2}=S.Fidelity;
	yLim{i,1}=nn(i)*Y_LIM;
	yLim{i,2}=Y_LIM;
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
			'XLim',X_LIM,'YLim',yLim{i,j},'XTick',X_TICK);
		fixticklength(TICK_LENGTH);
		if i~=M
			ax{i,j}.XTickLabel='';
		end
		if j~=1
			ax{i,j}.YAxisLocation='right';
		end
		for k=1:size(Y{i,j},1)
			plot(X{i},Y{i,j}(k,:),'Color',PLOT_LINE_COLOR{j}{k});
		end
	end
end

h=legend(ax{1,1},LEGEND_STRING{1},'NumColumns',2);
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
savefigure(fig,'F','qon','dynamics',2,P.CalcMode,...
	P.FileFormat,'Resolution',P.Resolution,'Rewrite',P.Rewrite);
close(fig);

end