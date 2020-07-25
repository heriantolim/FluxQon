function cID=SC_20180501_plot_exchange_time_vs_detuning_1(varargin)
%% Settings
Groot.usedefault();
Groot.usedefault('latex',12,.72);

PAPER_SIZE=[11.9,6.5];
AXES_POSITION=[1.5,1.3,6,5];

Y_LIM=[0,600];
TICK_LENGTH=.2;
X_TICK=-2:2:6;
X_TICK_LABEL=arrayfun(@(x)sprintf('%d$\\,\\omega$',x),...
	X_TICK,'UniformOutput',false);
X_TICK_LABEL{2}='0';
X_LABEL_STRING='Detuning, $\Delta$';
Y_LABEL_STRING='Exchange time (ns)';

LEGEND_STRING={'Large-$\Delta$ approx.','Calculated, n=1',...
	'Calculated, n=2','Calculated, n=3','Calculated, n=4'};

PLOT_LINE_COLOR='k';
PLOT_LINE_STYLE='--';

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
qb=runfunction({'default','qb'});
ion=runfunction({'default','ion'},1,'Concentration',P.Concentration);
mw=runfunction({'default','mw'});
gB=sqrt(ion.NumIons)*runfunction({'S','C',20170621,'get','Grms'},ion,qb);

F=FileDir('R','plot','exchange','time','vs','detuning',1,cID{:});
if exist([F.FullPath,'.mat'],'file')==2
	S1=loaddata(F);
	K=size(S1.Y,1);
else
	S1=struct();
	S2=loaddata('R','fidelity','vs','detuning',0,cID{:});
	K=numel(S2.NumPhotons);
	L=numel(S2.Detuning);
	S1.X=S2.Detuning;
	S1.Y=nan(K,L);
	for k=1:K
		l=1;
		[~,ix]=max(S2.Fidelity{k,l});
		S1.Y(k,l)=S2.Time{k,l}(ix);
		while l<L && S1.X(l)<0
			l=l+1;
			[Fx,ix]=max(S2.Fidelity{k,l}(S2.Time{k,l}<1.2*S1.Y(k,l-1)));
			if Fx<.99
				L1=l;
				break
			else
				S1.Y(k,l)=S2.Time{k,l}(ix);
			end
		end

		l=L;
		[~,ix]=max(S2.Fidelity{k,l});
		S1.Y(k,l)=S2.Time{k,l}(ix);
		while l>1 && S1.X(l)>0
			l=l-1;
			[Fx,ix]=max(S2.Fidelity{k,l}(S2.Time{k,l}<1.2*S1.Y(k,l+1)));
			if Fx<.99
				L2=l;
				break
			else
				S1.Y(k,l)=S2.Time{k,l}(ix);
			end
		end

		for l=L1:L2
			Fx=max(S2.Fidelity{k,l});
			if Fx>=.99
				[Fx,ix]=findpeaks(S2.Fidelity{k,l});
				m=find(Fx>=.99,1,'first');
				S1.Y(k,l)=S2.Time{k,l}(ix(m));
			end
		end
	end
	K=K+1;
	S1.Y(K,:)=NaN;
	ix=S1.X>=0;
	S1.Y(K,ix)=pi/2/gB^2*mw.Frequency*S1.X(ix);
	S1.Y=S1.Y*1e9;
	clear S2;
	savedata(S1,F,'Rewrite','yes');
end

xLim=S1.X([1,end]);

%% Drawing
fig=docfigure(PAPER_SIZE);

axes('Position',AXES_POSITION,...
	'XLim',xLim,'YLim',Y_LIM,'XTick',X_TICK,'XTickLabel',X_TICK_LABEL);
fixticklength(TICK_LENGTH);
xlabel(X_LABEL_STRING);
ylabel(Y_LABEL_STRING);

ix=41:50;
Y=S1.Y; Y(:,ix)=NaN;

p=cell(1,K);
for k=1:K-1
	p{k}=plot(S1.X,Y(k,:));
end
p{K}=plot(S1.X,S1.Y(K,:),'Color',PLOT_LINE_COLOR,'LineStyle',PLOT_LINE_STYLE);
uistack(p{K},'bottom');

for k=1:K-1
	l=[find(isnan(Y(k,:)),1,'first')-1,ix,find(isnan(Y(k,:)),1,'last')+1];
	plot(S1.X(l),S1.Y(k,l),'Marker','o','Color',p{k}.Color,'LineStyle','none');
end

h=legend(LEGEND_STRING,'Location','NorthEast');
p=h.Position;
p(1)=sum(AXES_POSITION([1,3]))+.2;
p(2)=sum(AXES_POSITION([2,4]))-p(4);
h.Position=p;

%% Saving
savefigure(fig,'F','exchange','time','vs','detuning','1',cID{:},...
	P.FileFormat,'Resolution',P.Resolution,'Rewrite',P.Rewrite);
close(fig);

end