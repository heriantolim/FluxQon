function SC_20160715_plot_state_and_frequency_1(varargin)
%% Settings
Groot.usedefault();
Groot.usedefault('latex',12,.72);

M=2; N=2;
L=Layout.tiledaxes(M,N,[5,5],[.1,.4,1.6,2],[0,.9,-.3,-.6]);
fprintf('Paper width  = %g\n',L.Paper.Position(3));
fprintf('Paper height = %g\n',L.Paper.Position(4));

X_LIM={[-pi/2,pi/2],[-pi/2,pi/2];[-2*pi,2*pi],[0,2*pi]};
Y_LIM={[-pi/2,pi/2],[-pi/2,pi/2]};
TICK_LENGTH=.2;
X_TICK={-pi/2:pi/4:pi/2,-pi/2:pi/4:pi/2;-2*pi:pi:2*pi,0:pi/2:2*pi};
X_TICK_LABEL={...
	{'-$\frac{\pi}{2}$','-$\frac{\pi}{4}$','0',...
		'$\frac{\pi}{4}$','$\frac{\pi}{2}$'},...
	{'-$\frac{\pi}{2}$','-$\frac{\pi}{4}$','0',...
		'$\frac{\pi}{4}$','$\frac{\pi}{2}$'};
	{'-2$\pi$','-$\pi$','0','$\pi$','2$\pi$'},...
	{'0','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','2$\pi$'}};
X_LABEL_STRING={'$\phi_+$','$\phi_+$';
	'Bias phase, $(\phi_\mathrm{bias}-\pi)\times 10^3$',...
	'Tuning phase, $\phi_\mathrm{tune}$'};
Y_LABEL_STRING={'$\phi_-$','$\phi_-$';
	'Excitation frequency, $\omega/\omega_0$',...
	'Excitation frequency, $\omega/\omega_0$'};
Y_LABEL_STRING2='Tunneling frequency, $\zeta/\omega_0$';
Z_LABEL_STRING={'Wavefunction, $\psi$','Wavefunction, $\psi$';'',''};
FIGURE_LABEL_STRING={'(a)','(b)';'(c)','(d)'};

LEGEND_STRING={'Excitation','Tunneling'};
LEGEND_LOCATION='East';

PLOT_LINE_COLOR={[220,80,20]/255,[0,120,190]/255};

%% Input Parsing
P=inputParser;
P.addParameter('Rewrite','yes',@(x)any(strcmpi(x,{'yes','no','ask'})));
P.addParameter('FileFormat',{'.fig','.pdf'},@FileDir.isfileext);
P.addParameter('Resolution',[300,600],@isintegervector);
P.parse(varargin{:});
P=P.Results;

%% Data Processing
W=loaddata('R','state','at','pi');
zLim1=[min([min(min(W(:,:,1))),min(min(W(:,:,2)))]),...
	max([max(max(W(:,:,1))),max(max(W(:,:,2)))])];
zLim1=zLim1*[1.25,-.25;-.25,1.25];
X=linspace(-pi/2,pi/2,size(W,1));

p=cell(1,N);
w=cell(1,N);
yLim2=cell(1,N);
cID={{'bias'},{'tuning'}};
for j=1:2
	S=loaddata('R','frequency','vs',cID{j}{:});
	p{j}=S(1,:);
	w{j}=S(2:3,:);
	yLim2{j}=[min(w{j},[],2),max(w{j},[],2)];
	w{j}(2,:)=diff(yLim2{j}(1,:))/diff(yLim2{j}(2,:)) ...
		*(w{j}(2,:)-yLim2{j}(2,1))+yLim2{j}(1,1);
end
p{1}=(p{1}-pi)*1e3;

%% Drawing
fig=docfigure(L.Paper.Position(3:4));

ax=cell(M,N);
for i=1:M
	for j=1:N
		ax{i,j}=axes('Position',L.Axes.Position{i,j},'XLim',X_LIM{i,j},...
			'XTick',X_TICK{i,j},'XTickLabel',X_TICK_LABEL{i,j});
		fixticklength(TICK_LENGTH);
		xlabel(X_LABEL_STRING{i,j});
		ylabel(Y_LABEL_STRING{i,j});
	end
end

for j=1:N
	set(ax{1,j},'ZLim',zLim1,'ZMinorTick','on',...
		'YLim',Y_LIM{1,j},'YTick',X_TICK{1,j},'YTickLabel',X_TICK_LABEL{1,j});
	view(ax{1,j},[1,5,1]);
	grid(ax{1,j},'on');
	zlabel(ax{1,j},Z_LABEL_STRING{1,j});
	set(ax{2,j},'Box','off','YLim',yLim2{j}(1,:));
end
set(ax{2,2},'YAxisLocation','right');

mesh(ax{1,2},X,X,W(:,:,2).');
ax{1,1}.CLim=ax{1,2}.CLim;
mesh(ax{1,1},X,X,W(:,:,1).');

for i=1:2
	for j=1:N
		plot(ax{2,j},p{j},w{j}(i,:),'Color',PLOT_LINE_COLOR{i});
	end
end

legend(ax{2,2},LEGEND_STRING,'Location',LEGEND_LOCATION);

oax=cell(1,N);
for j=1:2
	oax{j}=axes('Position',L.Axes.Position{2,j},...
		'Box','off','Color','none',...
		'XLim',X_LIM{2,j},'YLim',yLim2{j}(2,:),...
		'XTick',X_TICK{2,j},'XTickLabel',{},'XAxisLocation','top');
	fixticklength(TICK_LENGTH);
end
oax{1}.YAxisLocation='right';
ylabel(oax{1},Y_LABEL_STRING2);

%% Annotating
for i=1:M
	for j=1:N
		Label.subfigure(ax{i,j},FIGURE_LABEL_STRING{i,j},'Position',.4);
	end
end

%% Saving
savefigure(fig,'F','state','and','frequency',1,...
	P.FileFormat,'Resolution',P.Resolution,'Rewrite',P.Rewrite);
close(fig);

end
