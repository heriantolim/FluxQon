function SC_20180501_plot_g_effective_vs_tuning_flux_1(varargin)
%% Settings
Groot.usedefault();
Groot.usedefault('latex',12,.72);

PAPER_SIZE=[8.2,7.2];
AXES_POSITION=[1.5,1.3,5.5,5.5];

X_LIM=[.06,.17];
Y_LIM1=[0,600];
Y_LIM2=[0,6];
TICK_LENGTH=.2;
X_TICK=.05:.05:.2;
X_TICK_LABEL=arrayfun(@(x)sprintf('%.2f$\\,\\Phi_0$',x),...
	X_TICK,'UniformOutput',false);

Y_TICK1=0:100:600;

Y_TICK2=0:6;
Y_TICK_LABEL2=arrayfun(@(x)sprintf('%d$\\,\\omega_\\mathrm{mw}$',x),...
	Y_TICK2,'UniformOutput',false);
Y_TICK_LABEL2{1}='0';
Y_TICK_LABEL2{2}='$\omega_\mathrm{mw}$';

PLOT_LINE_COLOR1={[0,120,190]/255,[220,80,20]/255};
PLOT_LINE_STYLE1={'-','-'};
PLOT_LINE_COLOR2='k';
PLOT_LINE_STYLE2='-';


X_LABEL_STRING='Tuning flux, $\Phi_\mathrm{qb,tune}$';
Y_LABEL_STRING1='Coupling strength $/(2\pi)$ (kHz)';
Y_LABEL_STRING2='Qubit detuning';

LEGEND_STRING={...
	'$\mathrm{g}_{\mathrm{eff},15/2}$',...
	'$\mathrm{g}_{\mathrm{eff},13/2}$',...
	'$\omega_\mathrm{qb}-\omega_\mathrm{mw}$'};

%% Input Parsing
P=inputParser;
P.addParameter('Rewrite','yes',@(x)any(strcmpi(x,{'yes','no','ask'})));
P.addParameter('FileFormat',{'.fig','.pdf'},@FileDir.isfileext);
P.addParameter('Resolution',[300,600],@isintegervector);
P.parse(varargin{:});
P=P.Results;

%% Data Processing
qb=runfunction({'default','qb'});
ion=runfunction({'default','ion'},1);
ion.Multiplet=1:2;
G=runfunction({'S','C',20170621,'get','Grms'},ion,qb);

F=FileDir('R','plot','g','effective','vs','tuning','flux',1);
if exist([F.FullPath,'.mat'],'file')==2
	S1=loaddata(F);
else
	S1=struct();
	S2=loaddata('S','R',20160715,'frequency','vs','tuning');
	S3=loaddata('R','plot','osc','spectrum','vs','detuning',2,0,1,.005);
	
	ix=(S2(1,:)>=0 & S2(1,:)<pi) ...
		& (S2(2,:)-1>=S3.X{2}(1) & S2(2,:)-1<=S3.X{2}(end));
	S1.X=S2(1,ix);
	S1.Y1=[];
	S1.Y2=S2(2,ix)-1;
	S1.Y1=G*S2(3,ix)/S2(3,1).*interp1(S3.X{2},S3.Y1{2}(:,6).'/2,S1.Y2);
	
	% Clear memory.
	clear S2 S3;
	
	% Save data.
	savedata(S1,F,'Rewrite','yes');
end

S1.X=S1.X/2/pi;
S1.Y1=S1.Y1/2/pi;
S1.Y2=(S1.Y2-Y_LIM2(1))*diff(Y_LIM1)/diff(Y_LIM2)+Y_LIM1(1);

%% Drawing
fig=docfigure(PAPER_SIZE);

axes('Position',AXES_POSITION,'XLim',X_LIM,'YLim',Y_LIM2,...
	'TickLength',[0,0],'XTick',[],'YTick',Y_TICK2,...
	'YAxisLocation','right','YTickLabel',Y_TICK_LABEL2);
fixticklength(TICK_LENGTH);
ylabel(Y_LABEL_STRING2);

axes('Position',AXES_POSITION,'XLim',X_LIM,'YLim',Y_LIM1,...
	'XTick',X_TICK,'XTickLabel',X_TICK_LABEL,'YTick',Y_TICK1);
fixticklength(TICK_LENGTH);
xlabel(X_LABEL_STRING);
ylabel(Y_LABEL_STRING1);

for k=1:2
	plot(S1.X,S1.Y1(k,:),...
		'Color',PLOT_LINE_COLOR1{k},...
		'LineStyle',PLOT_LINE_STYLE1{k});
end

plot(S1.X,S1.Y2,...
	'Color',PLOT_LINE_COLOR2,...
	'LineStyle',PLOT_LINE_STYLE2);

legend(LEGEND_STRING,'Location','best');

%% Saving
savefigure(fig,'F','g','effective','vs','tuning','flux',1,...
	P.FileFormat,'Resolution',P.Resolution,'Rewrite',P.Rewrite);
close(fig);

end