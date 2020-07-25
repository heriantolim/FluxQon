function SC_20160715_plot_potential_1(O,varargin)
%% Settings
Groot.usedefault();
Groot.usedefault('latex',12,.72);

M=1; N=2;
L=Layout.tiledaxes(M,N,[5.5,5.5],[.1,.5,1,1.5],[.2,0,.3,0]);
fprintf('Paper width  = %g\n',L.Paper.Position(3));
fprintf('Paper height = %g\n',L.Paper.Position(4));

X_LIM={[-2*pi,2*pi],[-pi/2,pi/2]};
X_TICK={-2*pi:pi:2*pi,-pi/2:pi/4:pi/2};
Z_LIM=[-2,1];
TICK_LENGTH=.2;
X_TICK_LABEL={{'-2$\pi$','-$\pi$','0','$\pi$','2$\pi$'},...
	{'-$\frac{\pi}{2}$','-$\frac{\pi}{4}$','0',...
		'$\frac{\pi}{4}$','$\frac{\pi}{2}$'}};
X_LABEL_STRING={'Phase shift, $\phi_+$','$\phi_+$'};
Y_LABEL_STRING={'Phase shift, $\phi_-$','$\phi_-$'};
Z_LABEL_STRING='Potential, $U/E_J$';

FIGURE_LABEL_STRING={'(a)','(b)'};
FIGURE_LABEL_POSITION=[-1.3,0];

CONTOUR_LEVEL_STEP=.04;

%% Input Parsing
P=inputParser;
P.addRequired('NumPoints',@(x)isintegerscalar(x) && x>0);
P.addOptional('BiasPhase',pi,@isrealscalar);
P.addOptional('AlphaCoeff',.7,@(x)isrealscalar(x) && x>0);
P.addParameter('Rewrite','yes',@(x)any(strcmpi(x,{'yes','no','ask'})));
P.addParameter('FileFormat',{'.fig','.pdf'},@FileDir.isfileext);
P.addParameter('Resolution',[300,600],@isintegervector);
P.parse(O,varargin{:});
P=P.Results;

%% Data Processing
X=linspace(X_LIM{1}(1),X_LIM{1}(2),O);
U=zeros(O);
for k=1:O
	for l=1:O
		U(l,k)=-2*cos(X(k))*cos(X(l))-P.AlphaCoeff*cos(2*X(k)+P.BiasPhase);
	end
end
Umax=max(U(:));
Umin=min(U(:));

%% Drawing
fig=docfigure(L.Paper.Position(3:4));

ax=cell(M,N);
for i=1:M
	for j=1:N
		ax{i,j}=axes('Position',L.Axes.Position{i,j},...
			'XLim',X_LIM{i,j},'YLim',X_LIM{i,j},...
			'XTick',X_TICK{i,j},'YTick',X_TICK{i,j},...
			'XTickLabel',X_TICK_LABEL{i,j},'YTickLabel',X_TICK_LABEL{i,j});
		fixticklength(TICK_LENGTH);
		xlabel(X_LABEL_STRING{i,j});
		ylabel(Y_LABEL_STRING{i,j});
	end
end

contour(ax{1},X,X,U,'LevelStep',CONTOUR_LEVEL_STEP*(Umax-Umin));

ax{2}.CLim=ax{1}.CLim;
ax{2}.ZLim=Z_LIM;
zlabel(ax{2},Z_LABEL_STRING);
view(ax{2},[1,5,1]);
grid(ax{2},'on');

ix=X>=X_LIM{2}(1) & X<=X_LIM{2}(2);
mesh(ax{2},X(ix),X(ix),U(ix,ix));

%% Annotating
for i=1:M
	for j=1:N
		Label.subfigure(ax{i,j},FIGURE_LABEL_STRING{i,j},...
			'Position',FIGURE_LABEL_POSITION);
	end
end

%% Saving
savefigure(fig,'F','potential',1,...
	P.FileFormat,'Resolution',P.Resolution,'Rewrite',P.Rewrite);
close(fig);

end
