function cID=SC_20170626_compute_qon_mw_spectrum(varargin)

P=inputParser;
P.addOptional('Detuning',0,@isrealvector);% times Frequency
P.addOptional('IonBias',1,@isrealvector);% times B0
P.addOptional('CouplingRatio',1,@isrealscalar);
P.addOptional('Concentration',.005,@isrealscalar);% atomic percent
P.addOptional('ComponentID',{},@(x)isemptycell(x) || FileDir.iscomponentid(x));
P.addParameter('Frequency',2*pi*1.95e9,@isrealscalar);
P.addParameter('FockDimension',2,@isintegerscalar);
P.addParameter('TimeLength',1e-4,@isrealscalar);
P.addParameter('DispProgress',true);
P.parse(varargin{:});
P=P.Results;
cID=FileDir.expandlist(P.ComponentID);
cID=[{P.CouplingRatio,P.Concentration},cID];

qb=runfunction({'S','C',20180501,'default','qb'},'Frequency',P.Frequency);

ion=runfunction({'S','C',20180501,'default','ion'},1,...
	'Frequency',P.Frequency,'FockDimension',P.FockDimension,...
	'Concentration',P.Concentration);
a=round(ion.Rotation/pi*180);
G=runfunction({'S','C',20170621,'get','Grms'},...
	1,ion.Site,a(1),a(2),a(3),ion.Radius/qb.Radius,ion.Height/qb.Radius);

mw=runfunction({'S','C',20180501,'default','mw'},...
	'Frequency',P.Frequency,'FockDimension',P.FockDimension);
mw.MagneticAmplitude=P.CouplingRatio*G*sqrt(ion.NumIons) ...
	*Constant.FluxQuantum/pi/qb.Area/qb.Radius;

S=loaddata('S','R',20160715,'frequency','vs','tuning');
wz=1+P.Detuning;
wx=S(1,:)<pi;
wx=interp1(S(2,wx),S(3,wx),wz)*P.Frequency;
wz=wz*P.Frequency;
B0=ion.MagneticField;

S=struct();
S.Detuning=P.Detuning;
S.Frequency=[];
M=numel(wz);
for IB=P.IonBias
	ion.MagneticField=IB*B0;
	S.Intensity=cell(M,1);
	for i=1:M
		qb.Frequency=wz(i);
		qb.TunnelingFrequency=wx(i);
		ion.CouplingStrength=G*wx(i)/qb.Radius;
		RS=RadiationSpectrum(P.Frequency/pi,P.TimeLength,...
			Construct.Lindblad(qb,ion,mw),Hilbert.dimension(qb,ion,mw),3,...
			'DispProgress',P.DispProgress);
		S.Intensity{i}=RS(2,:);
	end
	S.Frequency=RS(1,:)*2*pi;
	S.Intensity=cell2mat(S.Intensity);
	savedata(S,'R','qon','mw','spectrum',IB,cID{:},'Rewrite','yes');
end

end