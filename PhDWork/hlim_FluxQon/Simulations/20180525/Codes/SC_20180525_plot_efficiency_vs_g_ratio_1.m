function SC_20180525_plot_efficiency_vs_g_ratio_1(varargin)
%% Settings
Groot.usedefault();
Groot.usedefault('latex',12,.72);

PAPER_SIZE=[12,5.8];
AXES_POSITION=[1.9,1.6,6,4];

Y_LIM=[0,1];
TICK_LENGTH=.2;
X_LABEL_STRING=['$\mathrm{log}_{10}\left(\frac',...
	'{\mathrm{g}_\mathrm{ion-mw}\sqrt{N}}',...
	'{\mathrm{g}_{\mathrm{ion-op}(u)}\sqrt{n}}\right)$'];
Y_LABEL_STRING=sprintf('Down-conversion\nefficiency');

%% Input Parsing
P=inputParser;
P.addOptional('CalcMode',0,@isintegerscalar);
P.addOptional('LineStrength3',1e3,@isrealscalar);
P.addOptional('ComponentID',{},@(x)isemptycell(x) || FileDir.iscomponentid(x));
P.addParameter('Rewrite','yes',@(x)any(strcmpi(x,{'yes','no','ask'})));
P.addParameter('FileFormat',{'.fig','.pdf'},@FileDir.isfileext);
P.addParameter('Resolution',[300,600],@isintegervector);
P.parse(varargin{:});
P=P.Results;
cID=FileDir.expandlist(P.ComponentID);

%% Data Processing
nn=[1,1,1,1,2,2,2];
NN=[1,2,3,1e4,2,3,1e4];
K=numel(NN);

X=cell(1,K);
Y=cell(1,K);
for k=1:K
	if NN(k)>10
		S2=loaddata('S','R',20180521,'ion','fidelity','vs','g','ratio',...
			P.CalcMode,nn(k),NN(k),P.LineStrength3,cID{:});
	else
		S2=loaddata('R','ion','fidelity','vs','g','ratio',...
			P.CalcMode,nn(k),NN(k),P.LineStrength3,cID{:});
	end
	L=numel(S2.CouplingRatio);
	X{k}=log10(sqrt(nn(k)/NN(k))*S2.CouplingRatio);
	Y{k}=nan(1,L);
	for l=1:L
		Y{k}(l)=max(S2.Fidelity{l});
		
		% Alternative method: efficiency = first maxima that is greater than .99
% 		Fx=max(S2.Fidelity{l});
% 		if Fx<.99
% 			Y{k}(l)=Fx;
% 		else
% 			Fx=findpeaks(S2.Fidelity{l});
% 			m=find(Fx>=.99,1,'first');
% 			Y{k}(l)=Fx(m);
% 		end
	end
end
clear S2;

xLim=[Inf,-Inf];
for k=1:K
	xLim(1)=min(xLim(1),X{k}(1));
	xLim(2)=max(xLim(2),X{k}(end));
end

legendString=cell(1,K);
for k=1:K
	if NN(k)>10
		legendString{k}=sprintf('$n=%d, N=10^{%d}$',nn(k),log10(NN(k)));
	else
		legendString{k}=sprintf('$n=%d, N=%d$',nn(k),NN(k));
	end
end

%% Drawing
fig=docfigure(PAPER_SIZE);

axes('Position',AXES_POSITION,'XLim',xLim,'YLim',Y_LIM);
xlabel(X_LABEL_STRING);
ylabel(Y_LABEL_STRING);
fixticklength(TICK_LENGTH);

for k=1:K
	plot(X{k},Y{k});
end

h=legend(legendString,'Location','NorthEast');
p=h.Position;
p(1)=sum(AXES_POSITION([1,3]))+.2;
p(2)=sum(AXES_POSITION([2,4]))-p(4);
h.Position=p;

%% Saving
savefigure(fig,'F','efficiency','vs','g','ratio',1,P.CalcMode,cID{:},...
	P.FileFormat,'Resolution',P.Resolution,'Rewrite',P.Rewrite);
close(fig);

end