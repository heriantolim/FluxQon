function SC_20180522_plot_ion_dynamics_1(varargin)
%% Settings
Groot.usedefault();
Groot.usedefault('latex',12,.72);

M=1; N=2;
L=Layout.tiledaxes(M,N,[6.4,4],[.1,.1,.1,.1],[.4,-.3,.7,.7],[.4,.4,.4,.4]);
fprintf('Paper width  = %g\n',L.Paper.Position(3));
fprintf('Paper height = %g\n',L.Paper.Position(4));

X_LIM={[0,500],[0,.1]};
Y_LIM=[-.1,1.3];
TICK_LENGTH=.2;
X_TICK={0:100:400,0:.02:.08};
X_LABEL_STRING='Time (ns)';
Y_LABEL_STRING='Excitation number';
FIGURE_LABEL_STRING={'(a) $N=10^8$','(b) $N=10^{13}$'};

LEGEND_STRING={...
	'$\langle\hat{b}_2^\dagger\hat{b}_2+\hat{b}_3^\dagger\hat{b}_3\rangle\quad$',...
	'$\langle\hat{c}_2^\dagger\hat{c}_2+\hat{c}_3^\dagger\hat{c}_3\rangle\quad$',...
	'$\langle\hat{a}_\mathrm{mw}^\dagger\hat{a}_\mathrm{mw}\rangle\quad$',....
	'$\langle\hat{a}_\mathrm{op}^\dagger\hat{a}_\mathrm{op}\rangle$'};

PLOT_LINE_COLOR={...
	[20,225,130]/255,...% green
	[20,190,230]/255,...% blue
	[125,50,150]/255,...% purple
	[255,130,40]/255};% orange

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
nn=1;
NN=[1e8,1e13];

X=cell(M,N);
Y=cell(M,N);
for i=1:M
	for j=1:N
		S=loaddata('R','ion','dynamics',P.CalcMode,nn,NN(j),.01,P.LineStrength3);
		X{i,j}=S.Time*1e9;
		Y{i,j}=[S.Number{2}+S.Number{3}; S.Number{5}+S.Number{6};
			S.Number{7}; S.Number{8}];
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
			'XLim',X_LIM{j},'YLim',Y_LIM*nn,'XTick',X_TICK{j});
		fixticklength(TICK_LENGTH);
		if i~=M
			ax{i,j}.XTickLabel='';
		end
		if j~=1
			ax{i,j}.YTickLabel='';
		end
		for k=1:size(Y{i,j},1)
			plot(X{i,j},Y{i,j}(k,:),'Color',PLOT_LINE_COLOR{k});
		end
	end
end
ax{1,1}.Children(2).LineStyle='--';

h=legend(ax{2},LEGEND_STRING,'NumColumns',4);
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
savefigure(fig,'F','ion','dynamics',1,P.CalcMode,...
	P.FileFormat,'Resolution',P.Resolution,'Rewrite',P.Rewrite);
close(fig);

end