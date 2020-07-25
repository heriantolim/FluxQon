function SC_20170621_plot_G_histogram_1(a,b,c,r,z,varargin)
%% Settings
Groot.usedefault();
Groot.usedefault('latex',12,.72);

PAPER_SIZE=[7.7,5.5];
AXES_POSITION=[1.4,1.3,6,4];

HIST_NUM_BINS=1e5;
PLOT_LINE_COLOR={[0,120,190]/255,[0,120,190]/255,...
	[220,80,20]/255,[220,80,20]/255};
PLOT_LINE_STYLE={'-','--','-','--'};
ANNOTATION_COLOR={'b','b','r','r'};
ANNOTATION_MARKER={'*','o','*','o'};

X_LIM=[0,10];
Y_LIM=[0,.2];
TICK_LENGTH=.2;
X_LABEL_STRING='g$_{\mathrm{qb-ion}}$ (kHz)';
Y_LABEL_STRING='Probability density';
LEGEND_LOCATION='NorthEast';
LEGEND_STRING={'J=15/2, Site 1','J=13/2, Site 1',...
	'J=15/2, Site 2','J=13/2, Site 2'};

%% Input Parsing
P=inputParser;
P.addRequired('Rot1Angle',@isrealscalar);
P.addRequired('Rot2Angle',@isrealscalar);
P.addRequired('Rot3Angle',@isrealscalar);
P.addRequired('CrystalRadius',@isrealscalar);% as a fraction of the qubit radius
P.addRequired('CrystalHeight',@isrealscalar);% as a fraction of the qubit radius
P.addOptional('QubitRadius',5e-6,@isrealscalar);
P.addOptional('TunnelingFrequency',6.375814463646857e+12,@isrealscalar);
P.addParameter('Rewrite','yes',@(x)any(strcmpi(x,{'yes','no','ask'})));
P.addParameter('FileFormat',{'.fig','.pdf'},@FileDir.isfileext);
P.addParameter('Resolution',[300,600],@isintegervector);
P.parse(a,b,c,r,z,varargin{:});
P=P.Results;

%% Data Processing
Hist_X=cell(2,2);
Hist_Y=cell(2,2);
r2=r^2;
S=loaddata('R','G',a,b,c);

for m=1:2
	for s=1:2
		S.G{m,s}=S.G{m,s}/P.QubitRadius*P.TunnelingFrequency/2/pi/1e3;
		Nx=numel(S.X);
		Ny=numel(S.Y);
		Nz=numel(S.Z);

		% Histogram
		Hist_Y{m,s}=nan(1,Nx*Ny*Nz);
		n=0;
		for i=1:Nx
			for j=1:Ny
				for k=1:Nz
					if S.X(i)^2+S.Y(j)^2<=r2 || S.Z(k)<=z
						n=n+1;
						Hist_Y{m,s}(n)=S.G{m,s}(i,j,k);
					end
				end
			end
		end
		Hist_Y{m,s}=Hist_Y{m,s}(1:n);
		[Hist_Y{m,s},Hist_X{m,s}]=histcounts(Hist_Y{m,s},HIST_NUM_BINS);
		Hist_X{m,s}=mean([Hist_X{m,s}(1:end-1);Hist_X{m,s}(2:end)],1);
		Hist_Y{m,s}=Hist_Y{m,s}/n;
		Hist_Y{m,s}(1)=0;
	end
end

%% Drawing
fig=docfigure(PAPER_SIZE);

axes('Position',AXES_POSITION,'XLim',X_LIM,'YLim',Y_LIM);
fixticklength(TICK_LENGTH);
xlabel(X_LABEL_STRING);
ylabel(Y_LABEL_STRING);

for k=1:4
	plot(Hist_X{k},Hist_Y{k},...
		'Color',PLOT_LINE_COLOR{k},...
		'LineStyle',PLOT_LINE_STYLE{k});
end

legend(LEGEND_STRING,'Location',LEGEND_LOCATION,'AutoUpdate','off');

%% Annotating
px=runfunction({'get','Grms'},1:2,1:2,a,b,c,r,z) ...
	/P.QubitRadius*P.TunnelingFrequency/2/pi/1e3;
py=zeros(size(px));
for k=1:4
	py(k)=mean(Hist_Y{k}(find(Hist_X{k}<px(k),1,'last')+[0,1]));
	plot(px(k),py(k),...
		'Color',ANNOTATION_COLOR{k},'Marker',ANNOTATION_MARKER{k});
end

[px,py]=data2pos(px,py);
annotation('textarrow','Position',[px(1),py(1)+.7,0,-.5],...
	'String','rms');

%% Saving
savefigure(fig,'F','G','histogram',1,a,b,c,r,z,...
	P.FileFormat,'Resolution',P.Resolution,'Rewrite',P.Rewrite);
close(fig);

end