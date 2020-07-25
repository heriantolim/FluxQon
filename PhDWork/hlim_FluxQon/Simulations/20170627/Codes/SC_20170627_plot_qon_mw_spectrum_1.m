function SC_20170627_plot_qon_mw_spectrum_1(varargin)
%% Settings
Groot.usedefault();
Groot.usedefault('latex',12,.72);

M=1; N=2;
L=Layout.tiledaxes(M,N,[6,6],[.1,.1,.5,.1],[-.5,-.7,.7,.7],[.5,.7,.1,.7]);
fprintf('Paper width  = %g\n',L.Paper.Position(3));
fprintf('Paper height = %g\n',L.Paper.Position(4));

X_LIM=[-3,3];
Y_LIM=[-4,8];
TICK_LENGTH=.2;
Y_TICK=[-5:2:-1,0,1:2:9];
Y_TICK_LABEL=arrayfun(@(x)sprintf('%d$\\,\\mathrm{g}_\\mathrm{B}$',x),...
	Y_TICK,'UniformOutput',false);
Y_TICK_LABEL(3:5)={'-$\mathrm{g}_\mathrm{B}$','0','$\mathrm{g}_\mathrm{B}$'};
X_LABEL_STRING='Qubit bias flux, $(\Phi_\mathrm{qb,bias}-\Phi_0/2)\times 10^4$';
Y_LABEL_STRING='Frequency, $\omega-\omega_\mathrm{A}$';
FIGURE_LABEL_STRING={'(a)','(b) $\mathrm{g}_\mathrm{A}=\mathrm{g}_\mathrm{B}$'};

MAX_INTENSITY=1e-5;
MAX_INTENSITY_COLOR=[53,0,211]/255;
PLOT_LINE_COLOR=[220,80,20]/255;
PLOT_LINE_STYLE='--';
ANNOTATION_COLOR=[45,175,255]/255;
ANNOTATION_LINE_STYLE='--';

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

F1=FileDir('R','plot','qon','mw','spectrum',1,P.CouplingRatio,P.Concentration);
if exist([F1.FullPath,'.mat'],'file')==2
	S1=loaddata(F1);
else
	S1=struct();
	S1.C1=cell(M,N);
	S1.X1=S1.C1;
	S1.Y1=S1.C1;
	S1.X2=S1.C1;
	S1.Y2=S1.C1;
	F2(1)=FileDir('R','ion','qb','spectrum',1,P.Concentration);
	F2(2)=FileDir('R','qon','mw','spectrum',1,P.CouplingRatio,P.Concentration);
	F2(3)=FileDir('S','R',20160715,'frequency','vs','bias');
	F2(4)=FileDir('R','qb','mw','spectrum',P.CouplingRatio,P.Concentration);
	c=rgb2hsv(MAX_INTENSITY_COLOR);
	for i=1:M
		for j=1:N
			S2=loaddata(F2(j));
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
		end
	end

	S2=loaddata(F2(3));
	S2(1,:)=(S2(1,:)/2/pi-.5)*1e4;
	S2(2,:)=(S2(2,:)-1)*P.Frequency/gB;
	ix=S2(1,:)>=X_LIM(1) & S2(1,:)<=X_LIM(2);
	S1.X2{1}=S2(1,ix);
	S1.Y2{1}=S2(2,ix);

	S2=loaddata(F2(4));
	S2.QubitBias=(S2.QubitBias-.5)*1e4;
	S2.Frequency=(S2.Frequency-P.Frequency)/gB;
	ix=find(S2.QubitBias>=X_LIM(1) & S2.QubitBias<=X_LIM(2));
	S1.X2{2}=S2.QubitBias(ix);
	iy=S2.Frequency<0;
	[~,l]=max(S2.Intensity(ix,iy),[],2);
	S1.Y2{2}(1,:)=S2.Frequency(l);
	for k=1:numel(S1.X2{2})
		iy=S2.Frequency>1.3*S1.X2{2}(k)^2;
		[~,l]=max(S2.Intensity(ix(k),iy),[],2);
		S1.Y2{2}(2,k)=S2.Frequency(l+find(iy,1,'first')-1);
	end

	clear S2;
	savedata(S1,F1,'Rewrite','yes');
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
[px,py]=data2pos(ax{1},[X_LIM(1),0],[-1,1]);
for j=1:2
	annotation('line','Position',[px(1),py(j),px(2)-px(1),0],...
		'Color',ANNOTATION_COLOR,'LineStyle',ANNOTATION_LINE_STYLE);
end

for i=1:M
	for j=1:N
		Label.subfigure(ax{i,j},FIGURE_LABEL_STRING{i,j},...
			'Location','North','Position',.4);
	end
end

%% Saving
savefigure(fig,'F','qon','mw','spectrum',1,P.CouplingRatio,P.Concentration,...
	P.FileFormat,'Resolution',P.Resolution,'Rewrite',P.Rewrite);
close(fig);

end