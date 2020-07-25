function SC_20170626_plot_qon_mw_spectrum_1(varargin)
%% Settings
R=[.5,1,2];

Groot.usedefault();
Groot.usedefault('latex',12,.72);

M=1; N=3;
L=Layout.tiledaxes(M,N,[4.5,4.5],[.1,.1,.5,.1],[-.4,-.7,.7,.7],[.5,.7,.1,.7]);
fprintf('Paper width  = %g\n',L.Paper.Position(3));
fprintf('Paper height = %g\n',L.Paper.Position(4));

X_LIM=[0,12];
Y_LIM=[-3,9];
TICK_LENGTH=.2;
X_TICK=0:3:12;
Y_TICK=-3:3:9;
X_TICK_LABEL=arrayfun(@(x)sprintf('%d$\\,\\mathrm{g}_\\mathrm{B}$',x),...
	X_TICK,'UniformOutput',false);
X_TICK_LABEL([1,end])={'0',''};
Y_TICK_LABEL=arrayfun(@(x)sprintf('%d$\\,\\mathrm{g}_\\mathrm{B}$',x),...
	Y_TICK,'UniformOutput',false);
Y_TICK_LABEL(2)={'0'};
X_LABEL_STRING='Qubit detuning, $\Delta$';
Y_LABEL_STRING='Frequency, $\omega-\omega_\mathrm{A}$';
FIGURE_LABEL_STRING={...
	'(a) $\mathrm{g}_\mathrm{A}=\frac{1}{2}\mathrm{g}_\mathrm{B}$',...
	'(b) $\mathrm{g}_\mathrm{A}=\mathrm{g}_\mathrm{B}$',...
	'(c) $\mathrm{g}_\mathrm{A}=2\mathrm{g}_\mathrm{B}$'};

MAX_INTENSITY=1e-5;
MAX_INTENSITY_COLOR=[53,0,211]/255;
PLOT_LINE_COLOR=[45,175,255]/255;
PLOT_LINE_STYLE='--';

%% Input Parsing
P=inputParser;
P.addOptional('IonBias',1,@isrealscalar);
P.addOptional('Concentration',.005,@isrealscalar);% atomic percent
P.addParameter('Frequency',2*pi*1.95e9,@isrealscalar);
P.addParameter('Rewrite','yes',@(x)any(strcmpi(x,{'yes','no','ask'})));
P.addParameter('FileFormat',{'.fig','.pdf'},@FileDir.isfileext);
P.addParameter('Resolution',[300,600],@isintegervector);
P.parse(varargin{:});
P=P.Results;

%% Data Processing
qb=runfunction({'S','C',20180501,'default','qb'},'Frequency',P.Frequency);
ion=runfunction({'S','C',20180501,'default','ion'},1,...
	'Frequency',P.Frequency,'Concentration',P.Concentration);
gB=runfunction({'S','C',20170621,'get','Grms'},ion,qb)*sqrt(ion.NumIons);

F=FileDir('R','plot','qon','mw','spectrum',1,P.IonBias,P.Concentration);
if exist([F.FullPath,'.mat'],'file')==2
	S1=loaddata(F);
else
	S1=struct();
	S1.C=cell(M,N);
	S1.X=S1.C;
	S1.Y=S1.C;
	c=rgb2hsv(MAX_INTENSITY_COLOR);
	for i=1:M
		for j=1:N
			S2=loaddata('R','qon','mw','spectrum',...
				P.IonBias,R(M*(i-1)+j),P.Concentration);
			S2.Detuning=S2.Detuning*P.Frequency/gB;
			S2.Frequency=(S2.Frequency-P.Frequency)/gB;
			ix=S2.Detuning>=X_LIM(1) & S2.Detuning<=X_LIM(2);
			iy=S2.Frequency>=Y_LIM(1) & S2.Frequency<=Y_LIM(2);
			S1.X{i,j}=downsample(S2.Detuning(ix),1000);
			S1.Y{i,j}=downsample(S2.Frequency(iy),1000);
			S1.C{i,j}=downsample(S2.Intensity(ix,iy).'/MAX_INTENSITY,1000);
			O=ones(size(S1.C{i,j}));
			S1.C{i,j}(S1.C{i,j}>1)=1;
			S1.C{i,j}=hsv2rgb(cat(3,c(1)*O,c(2)*S1.C{i,j},c(3)*O));
		end
	end

	clear S2;
	savedata(S1,F,'Rewrite','yes');
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
			'XLim',X_LIM,'YLim',Y_LIM,...
			'XTick',X_TICK,'YTick',Y_TICK,...
			'XTickLabel',X_TICK_LABEL,'YTickLabel',Y_TICK_LABEL);
		fixticklength(TICK_LENGTH);
		if i~=M
			ax{i,j}.XTickLabel='';
		end
		if j~=1
			ax{i,j}.YTickLabel='';
		end
		image(S1.X{i,j},S1.Y{i,j},S1.C{i,j});
		plot([0,Y_LIM(2)],[0,Y_LIM(2)],...
			'Color',PLOT_LINE_COLOR,'LineStyle',PLOT_LINE_STYLE);
	end
end

%% Annotating
for i=1:M
	for j=1:N
		Label.subfigure(ax{i,j},FIGURE_LABEL_STRING{i,j},'Position',.4);
	end
end

%% Saving
savefigure(fig,'F','qon','mw','spectrum',1,P.IonBias,P.Concentration,...
	P.FileFormat,'Resolution',P.Resolution,'Rewrite',P.Rewrite);
close(fig);

end