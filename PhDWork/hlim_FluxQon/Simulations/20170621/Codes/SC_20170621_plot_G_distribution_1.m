function SC_20170621_plot_G_distribution_1(m,a,b,c,r,z,varargin)
%% Settings
Groot.usedefault();
Groot.usedefault('latex',12,.72);

M=2; N=2;
L=Layout.tiledaxes(M,N,[4.5,4.5],[.8,.8,.8,.8],[-.5,1.1,.4,.3]);
fprintf('Paper width  = %g\n',L.Paper.Position(3));
fprintf('Paper height = %g\n',L.Paper.Position(4));

COLORBAR_SIZE=[.4,6];
CONTOUR_LEVEL_STEP=.1;
SHAPE_LINE_COLOR='r';

C_LIM=[0,12];
X_LIM=[-1,1];
Y_LIM={[0,2];[-1,1]};
TICK_LENGTH=.2;
X_TICK=-1:.5:1;
Y_TICK={0:.5:2;-1:.5:1};
X_TICK_LABEL={'-$R$','-$\frac{1}{2}R$','0','$\frac{1}{2}R$','$R$'};
Y_TICK_LABEL={...
	{'0','$\frac{1}{2}R$','$R$','$\frac{3}{2}R$','2$R$'};
	{'-$R$','-$\frac{1}{2}R$','0','$\frac{1}{2}R$','$R$'}};
C_LABEL_STRING='g$_{\mathrm{qb-ion},15/2}$ (kHz)';
X_LABEL_STRING='$x$';
Y_LABEL_STRING={'$z$';'$y$'};
Y_LABEL_ROTATION=0;
Y_LABEL_HORIZONTAL_ALIGNMENT={'right';'left'};
FIGURE_LABEL_STRING={'(a) Site 1','(b) Site 2';'(c)','(d)'};

%% Input Parsing
P=inputParser;
P.addRequired('Multiplet',@isintegerscalar);
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
P.parse(m,a,b,c,r,z,varargin{:});
P=P.Results;

%% Data Processing
Contour_X=cell(M,N);
Contour_Y=cell(M,N);
Contour_Z=cell(M,N);
r2=r^2;
S=loaddata('R','G',a,b,c);

for s=1:N
	S.G{m,s}=S.G{m,s}/P.QubitRadius*P.TunnelingFrequency/2/pi/1e3;
	Nx=numel(S.X);
	Ny=numel(S.Y);
	Nz=numel(S.Z);

	% Contour data along x-z plane
	Contour_X{1,s}=S.X;
	Contour_Y{1,s}=S.Z;
	j=find(S.Y>=0,1,'first');
	Contour_Z{1,s}=zeros(Nz,Nx);
	for i=1:Nx
		for k=1:Nz
			if S.X(i)<-r || S.X(i)>r || S.Z(k)>z
				Contour_Z{1,s}(k,i)=NaN;
			else
				Contour_Z{1,s}(k,i)=S.G{m,s}(i,j,k);
			end
		end
	end

	% Contour data along x-y plane
	Contour_X{2,s}=S.X;
	Contour_Y{2,s}=S.Y;
	k=find(S.Z>=.25*z,1,'first');
	Contour_Z{2,s}=zeros(Ny,Nx);
	for i=1:Nx
		for j=1:Ny
			if S.X(i)^2+S.Y(j)^2>r2
				Contour_Z{2,s}(j,i)=NaN;
			else
				Contour_Z{2,s}(j,i)=S.G{m,s}(i,j,k);
			end
		end
	end
end

%% Drawing
fig=docfigure(L.Paper.Position(3:4));

ax=cell(M,N);
for i=1:M
	for s=1:N
		ax{i,s}=axes('Position',L.Axes.Position{i,s},...
			'CLim',C_LIM,'XLim',X_LIM,'YLim',Y_LIM{i},...
			'XTick',X_TICK,'YTick',Y_TICK{i},...
			'XTickLabel',X_TICK_LABEL,'YTickLabel',Y_TICK_LABEL{i});
		fixticklength(TICK_LENGTH);
		xlabel(X_LABEL_STRING);
		ylabel(Y_LABEL_STRING{i},'Rotation',Y_LABEL_ROTATION,...
			'HorizontalAlignment',Y_LABEL_HORIZONTAL_ALIGNMENT{i});
	end
end

h=colorbar(ax{1,2},'Location','EastOutside');
ax{1,2}.Position=L.Axes.Position{1,2};
p=h.Position;
p(3:4)=COLORBAR_SIZE;
p(2)=mean([ax{1,2}.Position(2),sum(ax{2,2}.Position([2,4]))])-p(4)/2;
h.Position=p;
set(h.Label,'String',C_LABEL_STRING,'Interpreter','latex');

for i=1:M
	for s=1:N
		contour(ax{i,s},Contour_X{i,s},Contour_Y{i,s},Contour_Z{i,s},...
			'LevelStep',CONTOUR_LEVEL_STEP);
	end
end
for s=1:N
	rectangle(ax{1,s},'Position',[-r,0,2*r,z],...
		'EdgeColor',SHAPE_LINE_COLOR);
	rectangle(ax{2,s},'Position',[-r,-r,2*r,2*r],'Curvature',[1,1],...
		'EdgeColor',SHAPE_LINE_COLOR);
end

%% Annotating
for i=1:M
	for s=1:N
		Label.subfigure(ax{i,s},FIGURE_LABEL_STRING{i,s},'Position',.4);
	end
end

%% Saving
savefigure(fig,'F','G','distribution',1,m,a,b,c,r,z,...
	P.FileFormat,'Resolution',P.Resolution,'Rewrite',P.Rewrite);
close(fig);

end
