function SC_20180501_plot_osc_spectrum_vs_concentration_2(n,varargin)
%% Settings
Groot.usedefault();
Groot.usedefault('latex',12,.72);

PAPER_SIZE=[7,6.5];
AXES_POSITION=[1.9,1.3,5,5];

COLORBAR_SIZE=[.4,5];
ZERO_CUT_OFF=5e-5;% times Frequency
MIN_PEAK_DISTANCE=5;% times the delta f
MIN_PEAK_HEIGHT=.005;

X_LIM=[0,.05];
Y_LIM=[0,.04];
TICK_LENGTH=.2;
Y_TICK=0:.01:.04;
Y_TICK_LABEL=arrayfun(@(x)sprintf('%.2f$\\,\\omega$',x),...
	Y_TICK,'UniformOutput',false);
Y_TICK_LABEL{1}='0';

C_LABEL_STRING='Intensity';
X_LABEL_STRING='Ion concentration (at.\%)';
Y_LABEL_STRING='Frequency';

PLOT_LINE_COLOR='k';
PLOT_LINE_STYLE='--';

%% Input Parsing
P=inputParser;
P.addRequired('NumPhotons',@isintegerscalar);
P.addOptional('CalcMode',0,@isintegerscalar);
P.addOptional('Detuning',.5,@isrealscalar);% times Frequency
P.addOptional('CouplingRatio',1,@isrealscalar);
P.addParameter('Concentration',.005,@isrealscalar);% atomic percent
P.addParameter('Frequency',2*pi*1.95e9,@isrealscalar);
P.addParameter('Rewrite','yes',@(x)any(strcmpi(x,{'yes','no','ask'})));
P.addParameter('FileFormat',{'.fig','.pdf'},@FileDir.isfileext);
P.addParameter('Resolution',[300,600],@isintegervector);
P.parse(n,varargin{:});
P=P.Results;
cID={P.CalcMode,P.Detuning,P.CouplingRatio};

%% Data Processing
qb=runfunction({'default','qb'});
ion=runfunction({'default','ion'},1,'Concentration',P.Concentration);
gB=runfunction({'S','C',20170621,'get','Grms'},ion,qb)*sqrt(ion.NumIons);

F=FileDir('R','plot','osc','spectrum','vs','concentration',2,cID{:});
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

	% Load the data.
	S2=loaddata('R','osc','spectrum','vs','concentration',n,cID{:});
	S2.Frequency=S2.Frequency/P.Frequency;
	ix=S2.Concentration>=X_LIM(1) & S2.Concentration<=X_LIM(2);
	S2.Concentration=S2.Concentration(ix);
	S2.Intensity=S2.Intensity(ix,:);
	K=numel(S2.Concentration);

	if xlsWrite
		% Find peak locations.
		iy=S2.Frequency>ZERO_CUT_OFF;
		PL=zeros(K,1);
		for k=1:K
			Z=S2.Intensity(k,iy);
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
					lp=lp+1; PL(k,lp)=ix(l1+lx-1);
					l1=l3;
				end
				l2=l3;
			end
			[~,lx]=max(Z(ix(l1:l2)));
			lp=lp+1; PL(k,lp)=ix(l1+lx-1);
		end

		% Convert zeros to NaNs.
		PL(PL==0)=NaN;

		% Restore the starting index.
		PL=PL+(find(iy,1,'first')-1);

		% Save peak locations to the xls file.
		xlswrite(xlsPath,PL);

		% Wait for user to edit the xls file.
		input(sprintf(['Please edit %s\n',...
			'and align the connected frequency indices in columns,\n',...
			'with each column corresponding to a plot lineseries.\n',...
			'Then press Enter to continue ...'],[F.FileName,'.xls']),'s');
	end
	
	% Read the xls file.
	PL=xlsread(xlsPath);

	% Extract the data for plots.
	S1=struct();
	S1.X=S2.Concentration;
	S1.Y1=PL;
	S1.Z1=PL;
	ix=~isnan(PL);
	S1.Y1(ix)=S2.Frequency(PL(ix));
	for k=1:K
		S1.Z1(k,ix(k,:))=S2.Intensity(k,PL(k,ix(k,:)));
	end

	% 2 times large-detuning approximation of g_eff.
	ix=S1.X>0;
	S1.Y2=nan(1,K);
	S1.Y2(ix)=S1.X(ix)/P.Concentration/P.Detuning*2*gB^2/P.Frequency^2;

	% Clear memory.
	clear S2;
	
	% Save data.
	savedata(S1,F,'Rewrite','yes');
end

cLim=[0,max(S1.Z1(:))];

%% Drawing
fig=docfigure(PAPER_SIZE);

ax=axes('Position',AXES_POSITION,...
	'CLim',cLim,'Colormap',cool(256),...
	'XLim',X_LIM,'YLim',Y_LIM,'YTick',Y_TICK,'YTickLabel',Y_TICK_LABEL);
fixticklength(TICK_LENGTH);
xlabel(X_LABEL_STRING);
ylabel(Y_LABEL_STRING);

X=repmat(S1.X,2,1);
Z=repmat(zeros(size(S1.X)),2,1);
for k=1:size(S1.Y1,2)
	Y=repmat(S1.Y1(:,k).',2,1);
	C=repmat(S1.Z1(:,k).',2,1);
	surf(X,Y,Z,C,'FaceColor','none','EdgeColor','interp');
end

plot(S1.X,S1.Y2,'Color',PLOT_LINE_COLOR,'LineStyle',PLOT_LINE_STYLE);

% Color bar.
h=colorbar(ax,'Location','EastOutside');
ax.Position=AXES_POSITION;
p=h.Position;
p(3:4)=COLORBAR_SIZE;
p(2)=mean([AXES_POSITION(2),sum(AXES_POSITION([2,4]))])-p(4)/2;
h.Position=p;
set(h.Label,'String',C_LABEL_STRING,'Interpreter','latex');

%% Saving
savefigure(fig,'F','osc','spectrum','vs','concentration',2,cID{:},...
	P.FileFormat,'Resolution',P.Resolution,'Rewrite',P.Rewrite);
close(fig);

end