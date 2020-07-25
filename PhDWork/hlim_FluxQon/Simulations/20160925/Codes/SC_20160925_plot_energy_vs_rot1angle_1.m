function SC_20160925_plot_energy_vs_rot1angle_1(varargin)
%% Settings
Groot.usedefault();
Groot.usedefault('latex',12,.72);

M=1; N=2;
L=Layout.tiledaxes(M,N,[6,6],[.1,.1,.5,.5],[-.4,-.5,.7,.7],[.5,.5,.1,.1]);
fprintf('Paper width  = %g\n',L.Paper.Position(3));
fprintf('Paper height = %g\n',L.Paper.Position(4));

TICK_LENGTH=.2;
X_LABEL_STRING='$\gamma$ angle (degree)';
Y_LABEL_STRING='Energy (GHz)';
FIGURE_LABEL_STRING={'(a) Site 1','(b) Site 2'};

PLOT_LINE_COLOR={...
	{[228, 26, 28]/255,[ 55,126,184]/255,[ 77,175, 74]/255,[152, 78,163]/255,...
	 [255,127,  0]/255,[166, 86, 40]/255,[247,129,191]/255,[153,153,153]/255},...
	{'k'}};
PLOT_LINE_WIDTH=[.72,1];

%% Input Parsing
P=inputParser;
P.addOptional('MagneticField',.08,@isrealscalar);
P.addOptional('Rot2Angle',90,@isrealscalar);
P.addOptional('Rot3Angle',0,@isrealscalar);
P.addParameter('Class',1,@isintegerscalar);
P.addParameter('Rewrite','yes',@(x)any(strcmpi(x,{'yes','no','ask'})));
P.addParameter('FileFormat',{'.fig','.pdf'},@FileDir.isfileext);
P.addParameter('Resolution',[300,600],@isintegervector);
P.parse(varargin{:});
P=P.Results;

%% Data Processing
S=loaddata('R','energy','vs','rot1angle',...
	'in',P.MagneticField,'at',P.Rot2Angle,P.Rot3Angle);
X=S.Rot1Angle;
Y=S.Energy(:,:,P.Class);
for k=1:2
	for j=1:N
		Y{k,j}=Y{k,j}/Constant.Planck/1e9;
	end
end

%% Drawing
fig=docfigure(L.Paper.Position(3:4));

bgaxes('Position',L.Container.Position);
xlabel(X_LABEL_STRING);
ylabel(Y_LABEL_STRING);

ax=cell(M,N);
for i=1:M
	for j=1:N
		ax{i,j}=axes('Position',L.Axes.Position{i,j},'XLim',[X(1),X(end)]);
		fixticklength(TICK_LENGTH);
		if i~=M
			ax{i,j}.XTickLabel='';
		end
		if j~=1
			ax{i,j}.YTickLabel='';
		end
	end
end

yLim=zeros(N,2);
for j=1:N
	for k=1:2
		for l=1:size(Y{k,j},1)
			plot(ax{j},X,Y{k,j}(l,:),...
				'Color',PLOT_LINE_COLOR{k}{ceil(l/2)},...
				'LineWidth',PLOT_LINE_WIDTH(k));
		end
	end
	yLim(j,:)=ax{j}.YLim;
end
yLim=[min(yLim(:,1)),max(yLim(:,2))];
for j=1:N
	ax{j}.YLim=yLim;
end

for i=1:M
	for j=1:N
		Label.subfigure(ax{i,j},FIGURE_LABEL_STRING{i,j},'Position',.4);
	end
end

%% Saving
savefigure(fig,'F','energy','vs','rot1angle',1,...
	'in',P.MagneticField,'at',P.Rot2Angle,P.Rot3Angle,...
	P.FileFormat,'Resolution',P.Resolution,'Rewrite',P.Rewrite);
close(fig);

end
