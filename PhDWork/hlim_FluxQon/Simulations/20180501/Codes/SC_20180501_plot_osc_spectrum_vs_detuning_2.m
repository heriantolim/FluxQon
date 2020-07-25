function SC_20180501_plot_osc_spectrum_vs_detuning_2(n,varargin)
%% Settings
Groot.usedefault();
Groot.usedefault('latex',12,.72);

M=1; N=2;
LL=Layout.tiledaxes(M,N,[5,5],[.1,.1,.1,1.5],[-.3,1.6,.7,.8],[.4,.4,.4,-1]);
fprintf('Paper width  = %g\n',LL.Paper.Position(3));
fprintf('Paper height = %g\n',LL.Paper.Position(4));

COLORBAR_SIZE=[.4,5];
ZERO_CUT_OFF=0.01;% times g
MIN_PEAK_DISTANCE=5;% times the delta f
MIN_PEAK_HEIGHT=.005;

X_LIM={[-.25,.25],[.25,5]};
Y_LIM={[0,5],[0,.25]};
TICK_LENGTH=.2;
X_TICK={-.2:.2:.2,1:2:5};
Y_TICK={0:5,0:.05:.25};

X_TICK_LABEL{1}=arrayfun(@(x)sprintf('%.1f$\\,\\omega$',x),...
	X_TICK{1},'UniformOutput',false);
X_TICK_LABEL{1}{2}='0';
X_TICK_LABEL{2}=arrayfun(@(x)sprintf('%d$\\,\\omega$',x),...
	X_TICK{2},'UniformOutput',false);
X_TICK_LABEL{2}{1}='$\omega$';

Y_TICK_LABEL{1}=arrayfun(@(x)sprintf('%d$\\,\\mathrm{g}$',x),...
	Y_TICK{1},'UniformOutput',false);
Y_TICK_LABEL{1}{1}='0';
Y_TICK_LABEL{1}{2}='$\mathrm{g}$';
Y_TICK_LABEL{2}=arrayfun(@(x)sprintf('%.2f$\\,\\mathrm{g}$',x),...
	Y_TICK{2},'UniformOutput',false);
Y_TICK_LABEL{2}{1}='0';

C_LABEL_STRING='Intensity';
X_LABEL_STRING='Detuning, $\Delta$';
Y_LABEL_STRING='Frequency';
FIGURE_LABEL_STRING={'(a)','(b)';};

PLOT_LINE_COLOR='k';
PLOT_LINE_STYLE='--';

%% Input Parsing
P=inputParser;
P.addRequired('NumPhotons',@isintegerscalar);
P.addOptional('CalcMode',0,@isintegerscalar);
P.addOptional('CouplingRatio',1,@isrealscalar);
P.addOptional('Concentration',.005,@isrealscalar);% atomic percent
P.addParameter('Frequency',2*pi*1.95e9,@isrealscalar);
P.addParameter('Rewrite','yes',@(x)any(strcmpi(x,{'yes','no','ask'})));
P.addParameter('FileFormat',{'.fig','.pdf'},@FileDir.isfileext);
P.addParameter('Resolution',[300,600],@isintegervector);
P.parse(n,varargin{:});
P=P.Results;
cID={P.CalcMode,P.CouplingRatio,P.Concentration};

%% Data Processing
qb=runfunction({'default','qb'});
ion=runfunction({'default','ion'},1,'Concentration',P.Concentration);
gB=runfunction({'S','C',20170621,'get','Grms'},ion,qb)*sqrt(ion.NumIons);

F=FileDir('R','plot','osc','spectrum','vs','detuning',2,cID{:});
if exist([F.FullPath,'.mat'],'file')==2
	S1=loaddata(F);
else
	% xls file for edit.
	xlsPath=[F.FullPath,'.xls'];
	xlsWrite=true;

	if exist(xlsPath,'file')==2
		xlsWrite=strcmpi(input(sprintf(['%s alread exists.\n',...
			'Redo calculation and rewrite file? Y/N [N]:'],...
			[F.FileName,'.xls'])),'Y');
		if xlsWrite
			delete(xlsPath);
		end
	end
	
	PL=cell(1,N); S2=PL;
	K=zeros(1,N);

	% Load the data.
	for j=1:N
		S2{j}=loaddata('R','osc','spectrum','vs','detuning',n,cID{:},j);
		S2{j}.Frequency=S2{j}.Frequency/gB;
		ix=S2{j}.Detuning>=X_LIM{j}(1) & S2{j}.Detuning<=X_LIM{j}(2);
		S2{j}.Detuning=S2{j}.Detuning(ix);
		S2{j}.Intensity=S2{j}.Intensity(ix,:);
		K(j)=numel(S2{j}.Detuning);
	end

	if xlsWrite
		% Remake the xls file.
		for j=1:N
			% Find peak locations.
			iy=S2{j}.Frequency>ZERO_CUT_OFF;
			PL{j}=zeros(K(j),1);
			for k=1:K(j)
				Z=S2{j}.Intensity(k,iy);
				ix=find(Z>=MIN_PEAK_HEIGHT);
				L=numel(ix);
				if L==0
					continue
				end
				l1=1; l2=1; lp=0;
				while l2<L
					l3=l2+1;
					if ix(l3)-ix(l2)>MIN_PEAK_DISTANCE
						[~,lx]=max(Z(ix(l1:l2)));
						lp=lp+1; PL{j}(k,lp)=ix(l1+lx-1);
						l1=l3;
					end
					l2=l3;
				end
				[~,lx]=max(Z(ix(l1:l2)));
				lp=lp+1; PL{j}(k,lp)=ix(l1+lx-1);
			end

			% Convert zeros to NaNs.
			PL{j}(PL{j}==0)=NaN;

			% Restore the starting index.
			PL{j}=PL{j}+(find(iy,1,'first')-1);

			% Save peak locations to the xls file.
			xlswrite(xlsPath,PL{j},j);
		end

		% Wait for user to edit the xls file.
		input(sprintf(['Please edit %s\n',...
			'and align the connected frequency indices in columns,\n',...
			'with each column corresponding to a plot lineseries.\n',...
			'Then press Enter to continue ...'],[F.FileName,'.xls']),'s');
	end
	
	% Read the xls file.
	for j=1:N
		PL{j}=xlsread(xlsPath,j);
	end

	% Extract the data for plots.
	S1=struct();
	S1.X=cell(1,N); S1.Y1=S1.X; S1.Y2=S1.X; S1.Z1=S1.X;
	for j=1:N
		S1.X{j}=S2{j}.Detuning;
		S1.Y1{j}=PL{j};
		S1.Z1{j}=PL{j};
		ix=~isnan(PL{j});
		S1.Y1{j}(ix)=S2{j}.Frequency(PL{j}(ix));
		for k=1:K(j)
			S1.Z1{j}(k,ix(k,:))=S2{j}.Intensity(k,PL{j}(k,ix(k,:)));
		end
		
		% 2 times large-detuning approximation of g_eff.
		ix=S1.X{j}>0;
		S1.Y2{j}=nan(1,K(j));
		S1.Y2{j}(ix)=2*gB/P.Frequency./S1.X{j}(ix);
	end

	% Clear memory.
	clear S2;
	
	% Save data.
	savedata(S1,F,'Rewrite','yes');
end

cLim=0;
for j=1:N
	cLim=max(cLim,max(S1.Z1{j}(:)));
end
cLim=[0,cLim];

%% Drawing
fig=docfigure(LL.Paper.Position(3:4));

bgaxes('Position',LL.Container.Position);
xlabel(X_LABEL_STRING);

ax=cell(1,N);
for j=1:N
	ax{j}=axes('Position',LL.Axes.Position{j},...
		'CLim',cLim,'Colormap',cool(256),...
		'XLim',X_LIM{j},'YLim',Y_LIM{j},'XTick',X_TICK{j},'YTick',Y_TICK{j},...
		'XTickLabel',X_TICK_LABEL{j},'YTickLabel',Y_TICK_LABEL{j});
	fixticklength(TICK_LENGTH);
	
	X=repmat(S1.X{j},2,1);
	Z=repmat(zeros(size(S1.X{j})),2,1);
	for k=1:size(S1.Y1{j},2)
		Y=repmat(S1.Y1{j}(:,k).',2,1);
		C=repmat(S1.Z1{j}(:,k).',2,1);
		surf(X,Y,Z,C,'FaceColor','none','EdgeColor','interp');
	end
	
	plot(S1.X{j},S1.Y2{j},'Color',PLOT_LINE_COLOR,'LineStyle',PLOT_LINE_STYLE);
end

ylabel(ax{1},Y_LABEL_STRING);

% Color bar.
h=colorbar(ax{2},'Location','EastOutside');
ax{2}.Position=LL.Axes.Position{2};
p=h.Position;
p(3:4)=COLORBAR_SIZE;
p(2)=mean([ax{2}.Position(2),sum(ax{2}.Position([2,4]))])-p(4)/2;
h.Position=p;
set(h.Label,'String',C_LABEL_STRING,'Interpreter','latex');

%% Annotating
for i=1:M
	for j=1:N
		Label.subfigure(ax{i,j},FIGURE_LABEL_STRING{i,j},'Position',.4);
	end
end

%% Saving
savefigure(fig,'F','osc','spectrum','vs','detuning',2,cID{:},...
	P.FileFormat,'Resolution',P.Resolution,'Rewrite',P.Rewrite);
close(fig);

end