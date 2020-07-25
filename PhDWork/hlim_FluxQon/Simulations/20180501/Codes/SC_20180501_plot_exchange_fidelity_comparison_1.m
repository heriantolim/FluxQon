function SC_20180501_plot_exchange_fidelity_comparison_1(varargin)
%% Settings
Groot.usedefault();
Groot.usedefault('latex',12,.72);

PAPER_SIZE=[6.7,6.6];
AXES_POSITION=[1.4,1.3,5,5];

X_LIM=[0,5];
Y_LIM=[0,1];
TICK_LENGTH=.2;
X_TICK=1:4;
X_LABEL_STRING='Initial number of excitations, $n$';
Y_LABEL_STRING='Exchange fidelity';

LEGEND_STRING={'with the qubit','without the qubit'};
LEGEND_LOCATION='SouthWest';

PLOT_MARKER='o';
PLOT_MARKER_FACE_COLOR={'none',[220,80,20]/255};
PLOT_LINE_COLOR={[0,120,190]/255,[220,80,20]/255};

%% Input Parsing
P=inputParser;
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
N=1:4;
D=.5;
LS=[.06459,.06459,.06369,.06369];

Y=zeros(2,numel(N));
j=0;
for n=N
	j=j+1;
	S=loaddata('R','qon','mw','dynamics',n,2,D,1,P.Concentration);
	Y(1,j)=max(S.Fidelity);
	S=loaddata('R','ion','mw','dynamics',n,2,LS(j),P.Concentration);
	Y(2,j)=max(S.Fidelity);
end

%% Drawing
fig=docfigure(PAPER_SIZE);

axes('Position',AXES_POSITION,'XLim',X_LIM,'YLim',Y_LIM,'XTick',X_TICK);
xlabel(X_LABEL_STRING);
ylabel(Y_LABEL_STRING);
fixticklength(TICK_LENGTH);

p=cell(1,2);
for i=1:2
	p{i}=plot(N,Y(i,:),'Color',PLOT_LINE_COLOR{i},...
		'Marker',PLOT_MARKER,'MarkerFaceColor',PLOT_MARKER_FACE_COLOR{i});
end

legend(LEGEND_STRING,'Location',LEGEND_LOCATION,'AutoUpdate','off');

uistack(p{1},'up');

%% Saving
savefigure(fig,'F','exchange','fidelity','comparison',1,cID{:},...
	P.FileFormat,'Resolution',P.Resolution,'Rewrite',P.Rewrite);
close(fig);

end