function SC_20180501_plot_exchange_fidelity_vs_g_ratio_1(varargin)
%% Settings
Groot.usedefault();
Groot.usedefault('latex',12,.72);

M=1; N=2;
LL=Layout.tiledaxes(M,N,[6,4.5],[.1,.1,.1,.1],[.4,-.4,.7,1.1],[.3,.4,.3,.4]);
fprintf('Paper width  = %g\n',LL.Paper.Position(3));
fprintf('Paper height = %g\n',LL.Paper.Position(4));

Y_LIM=[0,1.25];
TICK_LENGTH=.2;
X_TICK=-.5:.5:1;
Y_TICK=0:.25:1;
X_LABEL_STRING='$\mathrm{log}_{10}(\mathrm{g}_\mathrm{A}/\mathrm{g}_\mathrm{B})$';
Y_LABEL_STRING='Exchange fidelity';
FIGURE_LABEL_STRING={'(a) without decoherence','(b) with decoherence'};
LEGEND_STRING={'$\Delta=0\quad$','$\Delta=\frac{1}{2}\omega\quad$',...
	'$\Delta=\omega$'};

%% Input Parsing
P=inputParser;
P.addOptional('CalcMode',0,@isintegerscalar);
P.addOptional('Concentration',.005,@isrealscalar);% atomic percent
P.addOptional('ComponentID',{},@(x)isemptycell(x) || FileDir.iscomponentid(x));
P.addParameter('Rewrite','yes',@(x)any(strcmpi(x,{'yes','no','ask'})));
P.addParameter('FileFormat',{'.fig','.pdf'},@FileDir.isfileext);
P.addParameter('Resolution',[300,600],@isintegervector);
P.parse(varargin{:});
P=P.Results;
cID=FileDir.expandlist(P.ComponentID);
cID=[{P.Concentration},cID];

%% Data Processing
D=[0,.5,1];
K=numel(D);
F=FileDir('R','plot','exchange','fidelity','vs','g','ratio',P.CalcMode,cID{:});
if exist([F.FullPath,'.mat'],'file')==2
	S1=loaddata(F);
else
	CM=P.CalcMode+[0,2];
	S1=struct();
	S1.X=cell(1,N);
	S1.Y=cell(M,N);
	for j=1:N
		S1.Y{1,j}=zeros(K,1);
		for k=1:K
			S2=loaddata('R','fidelity','vs','g','ratio',CM(j),D(k),cID{:});
			L=numel(S2.CouplingRatio);
			S1.Y{1,j}(k,L)=0;
			for l=1:L
				S1.Y{1,j}(k,l)=max(S2.Fidelity{1,l});
			end
		end
		S1.X{j}=log10(S2.CouplingRatio);
	end
	clear S2;
	savedata(S1,F,'Rewrite','yes');
end

xLim=S1.X{N}([1,end]);

%% Drawing
fig=docfigure(LL.Paper.Position(3:4));

bgaxes('Position',LL.Container.Position);
xlabel(X_LABEL_STRING);

ax=cell(M,N);
for i=1:M
	for j=1:N
		ax{i,j}=axes('Position',LL.Axes.Position{i,j},...
			'XLim',xLim,'YLim',Y_LIM,'XTick',X_TICK,'YTick',Y_TICK);
		if i~=M
			ax{i,j}.XTickLabel='';
		end
		if j~=1
			ax{i,j}.YTickLabel='';
		end
		fixticklength(TICK_LENGTH);
		for k=1:K
			plot(S1.X{j},S1.Y{i,j}(k,:));
		end
	end
end

ylabel(ax{1,1},Y_LABEL_STRING);

h=legend(ax{1,N},LEGEND_STRING,'NumColumns',K);
p=h.Position;
p(1)=mean([ax{1,1}.Position(1),sum(ax{1,2}.Position([1,3]))])-p(3)/2;
p(2)=max(sum(ax{1,1}.Position([2,4])),sum(ax{1,2}.Position([2,4])))+.2;
h.Position=p;

%% Annotating
for i=1:M
	for j=1:N
		Label.subfigure(ax{i,j},FIGURE_LABEL_STRING{i,j},'Position',.4);
	end
end

%% Saving
savefigure(fig,'F','exchange','fidelity','vs','g','ratio',1,P.CalcMode,cID{:},...
	P.FileFormat,'Resolution',P.Resolution,'Rewrite',P.Rewrite);
close(fig);

end