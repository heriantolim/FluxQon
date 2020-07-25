function SC_20180501_plot_exchange_fidelity_vs_detuning_2(varargin)
%% Settings
Groot.usedefault();
Groot.usedefault('latex',12,.72);

M=1; N=2;
LL=Layout.tiledaxes(M,N,[6,5],[.1,.1,.1,1.7],[.4,-.4,.7,.9],[.4,.4,.4,-1.2]);
fprintf('Paper width  = %g\n',LL.Paper.Position(3));
fprintf('Paper height = %g\n',LL.Paper.Position(4));

Y_LIM={[0,1.05],[0,600]};
TICK_LENGTH=.2;
X_TICK=-2:2:6;
X_TICK_LABEL=arrayfun(@(x)sprintf('%d$\\,\\omega$',x),...
	X_TICK,'UniformOutput',false);
X_TICK_LABEL{2}='0';
X_LABEL_STRING='Detuning, $\Delta$';
Y_LABEL_STRING={'Exchange fidelity';'Exchange time (ns)'};
FIGURE_LABEL_STRING={'(a)','(b)'};
FIGURE_LABEL_LOCATION={'NorthEast','NorthWest'};

%% Input Parsing
P=inputParser;
P.addOptional('CouplingRatio',1,@isrealscalar);
P.addOptional('Concentration',.005,@isrealscalar);% atomic percent
P.addOptional('ComponentID',{},@(x)isemptycell(x) || FileDir.iscomponentid(x));
P.addParameter('Rewrite','yes',@(x)any(strcmpi(x,{'yes','no','ask'})));
P.addParameter('FileFormat',{'.fig','.pdf'},@FileDir.isfileext);
P.addParameter('Resolution',[300,600],@isintegervector);
P.parse(varargin{:});
P=P.Results;
cID=FileDir.expandlist(P.ComponentID);
cID=[{P.CouplingRatio,P.Concentration},cID];

%% Data Processing
F=FileDir('R','plot','exchange','fidelity','vs','detuning',2,cID{:});
if exist([F.FullPath,'.mat'],'file')==2
	S1=loaddata(F);
	K=numel(S1.N);
else
	S1=struct();
	S2=loaddata('R','fidelity','vs','detuning',2,cID{:});
	S1.N=S2.NumPhotons;
	S1.X=S2.Detuning;
	K=numel(S1.N);
	L=numel(S1.X);
	S1.Y=repmat({nan(K,L)},1,2);
	for k=1:K
		l=1;
		[Fx,ix]=max(S2.Fidelity{k,l});
		S1.Y{1}(k,l)=Fx;
		S1.Y{2}(k,l)=S2.Time{k,l}(ix);
		while l<L && S1.X(l)<0
			l=l+1;
			[Fx,ix]=max(S2.Fidelity{k,l}(S2.Time{k,l}<1.2*S1.Y{2}(k,l-1)));
			if Fx<.99
				L1=l;
				break
			else
				S1.Y{1}(k,l)=Fx;
				S1.Y{2}(k,l)=S2.Time{k,l}(ix);
			end
		end

		l=L;
		[Fx,ix]=max(S2.Fidelity{k,l});
		S1.Y{1}(k,l)=Fx;
		S1.Y{2}(k,l)=S2.Time{k,l}(ix);
		while l>1 && S1.X(l)>0
			l=l-1;
			[Fx,ix]=max(S2.Fidelity{k,l}(S2.Time{k,l}<1.2*S1.Y{2}(k,l+1)));
			if Fx<.99
				L2=l;
				break
			else
				S1.Y{1}(k,l)=Fx;
				S1.Y{2}(k,l)=S2.Time{k,l}(ix);
			end
		end

		for l=L1:L2
			[S1.Y{1}(k,l),ix]=max(S2.Fidelity{k,l});
			S1.Y{2}(k,l)=S2.Time{k,l}(ix);
		end
	end
	S1.Y{2}=S1.Y{2}*1e9;
	clear S2;
	savedata(S1,F,'Rewrite','yes');
end

xLim=S1.X([1,end]);

legendString=cell(1,K);
for k=1:K
	if k<K
		legendString{k}=sprintf('$n=%d\\quad$',S1.N(k));
	else
		legendString{k}=sprintf('$n=%d$',S1.N(k));
	end
end

%% Drawing
fig=docfigure(LL.Paper.Position(3:4));

bgaxes('Position',LL.Container.Position);
xlabel(X_LABEL_STRING);

ax=cell(M,N);
for i=1:M
	for j=1:N
		ax{i,j}=axes('Position',LL.Axes.Position{i,j},...
			'XLim',xLim,'YLim',Y_LIM{i,j},...
			'XTick',X_TICK,'XTickLabel',X_TICK_LABEL);
		fixticklength(TICK_LENGTH);
		for k=1:K
			plot(S1.X,S1.Y{i,j}(k,:));
		end
	end
end

ylabel(ax{1},Y_LABEL_STRING{1});
ylabel(ax{2},Y_LABEL_STRING{2});

h=legend(ax{1},legendString,'NumColumns',K);
p=h.Position;
p(1)=mean([ax{1}.Position(1),sum(ax{2}.Position([1,3]))])-p(3)/2;
p(2)=max(sum(ax{1}.Position([2,4])),sum(ax{2}.Position([2,4])))+.3;
h.Position=p;

%% Annotating
for i=1:M
	for j=1:N
		Label.subfigure(ax{i,j},FIGURE_LABEL_STRING{i,j},...
			'Position',.4,'Location',FIGURE_LABEL_LOCATION{i,j});
	end
end

%% Saving
savefigure(fig,'F','exchange','fidelity','vs','detuning',2,cID{:},...
	P.FileFormat,'Resolution',P.Resolution,'Rewrite',P.Rewrite);
close(fig);

end