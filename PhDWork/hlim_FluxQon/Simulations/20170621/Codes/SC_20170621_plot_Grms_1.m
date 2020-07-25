function SC_20170621_plot_Grms_1(r,z,varargin)
%% Settings
Groot.usedefault();
Groot.usedefault('latex',12,.72);
Groot.setdefault(...
	'shapeColor$',[70,70,70]/255,...
	'shapeLineStyle$',':');

M=2; N=2;
L=Layout.tiledaxes(M,N,[7,4],[.1,.1,.4,.1],[-.15,-.15,.7,1.05],[.4,.1,.1,.1]);
fprintf('Paper width  = %g\n',L.Paper.Position(3));
fprintf('Paper height = %g\n',L.Paper.Position(4));

X_LIM=[0,180];
Y_LIM={[0,25],[0,8]};
TICK_LENGTH=.2;
X_LABEL_STRING='$\gamma$ angle (degree)';
Y_LABEL_STRING={'Frequency $/(2\pi)$ (GHz)';
	'g$^\mathrm{rms}_\mathrm{qb-ion}/(2\pi)$ (kHz)'};
FIGURE_LABEL_STRING={...
	'(a) ($\gamma,0^\circ,0^\circ$)',...
	'(b) ($\gamma,-90^\circ,0^\circ$)';'(c)','(d)'};

LEGEND_LOCATION={'NorthEast';'NorthEast'};
LEGEND_STRING={...
	{'$\nu_1-\nu_0$ (Site 1)','$\nu_3-\nu_2$ (Site 1)',...
	'$\nu_1-\nu_0$ (Site 2)','$\nu_3-\nu_2$ (Site 2)'};
	{'$J=15/2$, Site 1','$J=13/2$, Site 1',...
	'$J=15/2$, Site 2','$J=13/2$, Site 2'}};

PLOT_LINE_COLOR={...
	{[255,125,0]/255,[255,125,0]/255,[30,180,30]/255,[30,180,30]/255};
	{[0,120,190]/255,[0,120,190]/255,[220,80,20]/255,[220,80,20]/255}};
PLOT_LINE_STYLE={'-','--','-','--'};

%% Input Parsing
P=inputParser;
P.addRequired('CrystalRadius',@isrealscalar);% as a fraction of the qubit radius
P.addRequired('CrystalHeight',@isrealscalar);% as a fraction of the qubit radius
P.addOptional('QubitRadius',5e-6,@isrealscalar);
P.addOptional('TunnelingFrequency',6.375814463646857e+12,@isrealscalar);
P.addParameter('ResonanceFrequency',2*pi*1.95e9,@isrealscalar);
P.addParameter('Rewrite','yes',@(x)any(strcmpi(x,{'yes','no','ask'})));
P.addParameter('FileFormat',{'.fig','.pdf'},@FileDir.isfileext);
P.addParameter('Resolution',[300,600],@isintegervector);
P.parse(r,z,varargin{:});
P=P.Results;
fR=P.ResonanceFrequency/2/pi/1e9;

%% Data Processing
b=[0,270];
A=0:180;
E=cell(1,N);
G=cell(1,N);
for j=1:N
	E{j}=runfunction({'get','frequency'},1:2,1:2,A,b(j),0)/2/pi/1e9;
	G{j}=runfunction({'get','Grms'},1:2,1:2,A,b(j),0,r,z) ...
		/P.QubitRadius*P.TunnelingFrequency/2/pi/1e3;
end

a=zeros(1,2);
gM=a;
for s=1:2
	[gM(s),a(s)]=max(G{2}(1,s,:));
end
[gM,s]=max(gM);
B=fR/E{2}(1,s,a(s));
a=A(a(s));

for j=1:N
	E{j}=reshape(E{j},4,[])*B;
	G{j}=reshape(G{j},4,[]);
end

%% Drawing
fig=docfigure(L.Paper.Position(3:4));

bgaxes('Position',L.Container.Position);
xlabel(X_LABEL_STRING);

ax=cell(M,N);
for i=1:M
	for j=1:N
		ax{i,j}=axes('Position',L.Axes.Position{i,j},...
			'XLim',X_LIM,'YLim',Y_LIM{i});
		fixticklength(TICK_LENGTH);
		if i~=M
			ax{i,j}.XTickLabel='';
		end
		if j~=1
			ax{i,j}.YTickLabel='';
		end
	end
	ylabel(ax{i,1},Y_LABEL_STRING{i});
end

for j=1:N
	for k=1:4
		plot(ax{1,j},A,E{j}(k,:),...
			'Color',PLOT_LINE_COLOR{1}{k},...
			'LineStyle',PLOT_LINE_STYLE{k});
		plot(ax{2,j},A,G{j}(k,:),...
			'Color',PLOT_LINE_COLOR{2}{k},...
			'LineStyle',PLOT_LINE_STYLE{k});
	end
end

for i=1:M
	legend(ax{i,1},LEGEND_STRING{i},'Location',LEGEND_LOCATION{i});
end

%% Annotating
px=cell(2,1);
py=cell(2,1);
[px{1},py{1}]=data2pos(ax{1,2},[0,a],[fR,0]);
[px{2},py{2}]=data2pos(ax{2,2},[0,a],[gM,ax{2,2}.YLim(2)]);

for i=1:M
	annotation('line','Position',[px{i}(1),py{i}(1),px{i}(2)-px{i}(1),0]);
	annotation('line','Position',[px{i}(2),py{i}(2),0,py{i}(1)-py{i}(2)]);
end

p=zeros(1,4);
p(3)=1;
p(4)=.4;
p(1)=mean([px{1}(2),px{2}(2)])-.5*p(3);
p(2)=mean([ax{1,2}.Position(2),sum(ax{2,2}.Position([2,4]))])-.5*p(4)-.06;
annotation('textbox','Position',p,'String',sprintf('%d',a));

p(1)=ax{1,2}.Position(1)+.2;
p(2)=py{1}(1)-.06;
annotation('textbox','Position',p,'String',sprintf('%.3f',fR));

p(1)=ax{2,2}.Position(1)+.2;
p(2)=py{2}(1)-p(4)-.07;
annotation('textbox','Position',p,'String',sprintf('%.3f',gM));

for i=1:M
	for j=1:N
		Label.subfigure(ax{i,j},FIGURE_LABEL_STRING{i,j},'Position',.4);
	end
end

%% Saving
savefigure(fig,'F','Grms',1,r,z,...
	P.FileFormat,'Resolution',P.Resolution,'Rewrite',P.Rewrite);
close(fig);

end
