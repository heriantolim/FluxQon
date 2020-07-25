function SC_20170627_plot_qon_mw_spectrum_2(varargin)
%% Settings
IB=[.91,.94,.97,1.03,1.06,1.09];

Groot.usedefault();
Groot.usedefault('latex',12,.72);

M=2; N=3;
L=Layout.tiledaxes(M,N,[4.5,4.5],[.1,.1,.5,.1],[-.5,-.7,.7,.7],[.5,.7,.1,.7]);
fprintf('Paper width  = %g\n',L.Paper.Position(3));
fprintf('Paper height = %g\n',L.Paper.Position(4));

X_LIM=[-3,3];
Y_LIM=[-4,8];
TICK_LENGTH=.2;
Y_TICK=-5:2:9;
Y_TICK_LABEL=arrayfun(@(x)sprintf('%d$\\,\\mathrm{g}_\\mathrm{B}$',x),...
	Y_TICK,'UniformOutput',false);
Y_TICK_LABEL(3:4)={'-$\mathrm{g}_\mathrm{B}$','$\mathrm{g}_\mathrm{B}$'};
X_LABEL_STRING='Qubit bias flux, $(\Phi_\mathrm{qb,bias}-\Phi_0/2)\times 10^4$';
Y_LABEL_STRING='Frequency, $\omega-\omega_\mathrm{A}$';
FIGURE_LABEL_STRING=arrayfun(@(k)sprintf(...
	'(%s) $B_\\mathrm{ion,bias}=%g\\,B_\\mathrm{R}$',96+k,IB(k)),...
	1:numel(IB),'UniformOutput',false);

MAX_INTENSITY=1e-5;
MAX_INTENSITY_COLOR=[53,0,211]/255;
PLOT_LINE_COLOR=[30,160,20]/255;
PLOT_LINE_STYLE='--';

%% Input Parsing
P=inputParser;
P.addOptional('CouplingRatio',1,@isrealscalar);
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

F=FileDir('R','plot','qon','mw','spectrum',2,P.CouplingRatio,P.Concentration);
if exist([F.FullPath,'.mat'],'file')==2
	S1=loaddata(F);
else
	S1=struct();
	S1.C1=cell(M,N);
	S1.X1=S1.C1;
	S1.Y1=S1.C1;
	S1.X2=S1.C1;
	S1.Y2=S1.C1;
	c=rgb2hsv(MAX_INTENSITY_COLOR);
	for i=1:M
		for j=1:N
			k=N*(i-1)+j;

			S2=loaddata('R','qon','mw','spectrum',...
				IB(k),P.CouplingRatio,P.Concentration);
			S2.QubitBias=(S2.QubitBias-.5)*1e4;
			S2.Frequency=(S2.Frequency-P.Frequency)/gB;
			ix=S2.QubitBias>=X_LIM(1) & S2.QubitBias<=X_LIM(2);
			iy=S2.Frequency>=Y_LIM(1) & S2.Frequency<=Y_LIM(2);
			S1.X1{i,j}=downsample(S2.QubitBias(ix),1000);
			S1.Y1{i,j}=downsample(S2.Frequency(iy),1000);
			S1.C1{i,j}=downsample(S2.Intensity(ix,iy).'/MAX_INTENSITY,1000);
			O=ones(size(S1.C1{i,j}));
			S1.C1{i,j}(S1.C1{i,j}>1)=1;
			S1.C1{i,j}=hsv2rgb(cat(3,c(1)*O,c(2)*S1.C1{i,j},c(3)*O));

			ix=S2.QubitBias>-2 & S2.QubitBias<2;
			switch k
				case {1,2}
					iy=S2.Frequency<(IB(k)-1)*P.Frequency/gB;
					[~,l]=max(S2.Intensity(ix,iy),[],2);
					S1.X2{i,j}=S2.QubitBias(ix);
					S1.Y2{i,j}=S2.Frequency(l);
				case {5,6}
					iy=S2.Frequency>(IB(k)-1)*P.Frequency/gB;
					[~,l]=max(S2.Intensity(ix,iy),[],2);
					S1.X2{i,j}=S2.QubitBias(ix);
					S1.Y2{i,j}=S2.Frequency(l+find(iy,1,'first')-1);
			end
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
			'XLim',X_LIM,'YLim',Y_LIM,'YTick',Y_TICK,'YTickLabel',Y_TICK_LABEL);
		fixticklength(TICK_LENGTH);
		if i~=M
			ax{i,j}.XTickLabel='';
		end
		if j~=1
			ax{i,j}.YTickLabel='';
		end
		image(S1.X1{i,j},S1.Y1{i,j},S1.C1{i,j});
		plot(S1.X2{i,j},S1.Y2{i,j},...
			'Color',PLOT_LINE_COLOR,'LineStyle',PLOT_LINE_STYLE);
	end
end

%% Annotating
for i=1:M
	for j=1:N
		Label.subfigure(ax{i,j},FIGURE_LABEL_STRING{N*(i-1)+j},...
			'Location','North','Position',.4);
	end
end

%% Saving
savefigure(fig,'F','qon','mw','spectrum',2,P.CouplingRatio,P.Concentration,...
	P.FileFormat,'Resolution',P.Resolution,'Rewrite',P.Rewrite);
close(fig);

end