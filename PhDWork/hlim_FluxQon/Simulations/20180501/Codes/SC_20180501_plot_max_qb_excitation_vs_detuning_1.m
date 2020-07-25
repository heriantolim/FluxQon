function SC_20180501_plot_max_qb_excitation_vs_detuning_1(varargin)
%% Settings
Groot.usedefault();
Groot.usedefault('latex',12,.72);

PAPER_SIZE=[8,6];
AXES_POSITION=[1.9,1.3,6,4.5];

Y_LIM=[0,1];
TICK_LENGTH=.2;
X_TICK=-2:2:6;
X_TICK_LABEL=arrayfun(@(x)sprintf('%d$\\,\\omega$',x),...
	X_TICK,'UniformOutput',false);
X_TICK_LABEL{2}='0';
X_LABEL_STRING='Detuning, $\Delta$';
Y_LABEL_STRING=sprintf('Qubit oscillation\namplitude');

LEGEND_LOCATION='NorthEast';

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
F=FileDir('R','plot','max','qb','excitation','vs','detuning',1,cID{:});
if exist([F.FullPath,'.mat'],'file')==2
	S1=loaddata(F);
	K=numel(S1.N);
else
	S1=struct();
	S2=loaddata('R','fidelity','vs','detuning',0,cID{:});
	S1.N=S2.NumPhotons;
	S1.X=S2.Detuning;
	K=numel(S1.N);
	L=numel(S1.X);
	S1.Y=nan(K,L);
	for k=1:K
		for l=1:L
			S1.Y(k,l)=max(S2.Number{k,l});
		end
	end
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
fig=docfigure(PAPER_SIZE);

axes('Position',AXES_POSITION,'XLim',xLim,'YLim',Y_LIM,...
	'XTick',X_TICK,'XTickLabel',X_TICK_LABEL);
xlabel(X_LABEL_STRING);
ylabel(Y_LABEL_STRING);
fixticklength(TICK_LENGTH);

for k=1:K
	plot(S1.X,S1.Y(k,:));
end

legend(legendString,'Location',LEGEND_LOCATION);

%% Saving
savefigure(fig,'F','max','qb','excitation','vs','detuning',1,cID{:},...
	P.FileFormat,'Resolution',P.Resolution,'Rewrite',P.Rewrite);
close(fig);

end