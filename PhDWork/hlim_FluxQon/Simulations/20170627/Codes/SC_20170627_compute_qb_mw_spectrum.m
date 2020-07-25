function cID=SC_20170627_compute_qb_mw_spectrum(varargin)

P=inputParser;
P.addOptional('QubitBias',.5,@isrealvector);% times FluxQuantum
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
	'Frequency',P.Frequency,'Concentration',P.Concentration);
a=round(ion.Rotation/pi*180);
G=runfunction({'S','C',20170621,'get','Grms'},...
	1,ion.Site,a(1),a(2),a(3),ion.Radius/qb.Radius,ion.Height/qb.Radius);

mw=runfunction({'S','C',20180501,'default','mw'},...
	'Frequency',P.Frequency,'FockDimension',P.FockDimension);
mw.MagneticAmplitude=P.CouplingRatio*G*sqrt(ion.NumIons) ...
	*Constant.FluxQuantum/pi/qb.Area/qb.Radius;

w=loaddata('S','R',20160715,'frequency','vs','bias').';
w=interp1(w(:,1),w(:,2:3),2*pi*(P.QubitBias-floor(P.QubitBias)))*P.Frequency;

S=struct();
S.QubitBias=P.QubitBias;
S.Frequency=[];
M=size(w,1);
S.Intensity=cell(M,1);
for i=1:M
	qb.Frequency=w(i,1);
	qb.TunnelingFrequency=w(i,2);
	RS=RadiationSpectrum(P.Frequency/pi,P.TimeLength,...
		Construct.Lindblad(qb,mw),Hilbert.dimension(qb,mw),2,...
		'DispProgress',P.DispProgress);
	S.Intensity{i}=RS(2,:);
end
S.Frequency=RS(1,:)*2*pi;
S.Intensity=cell2mat(S.Intensity);
savedata(S,'R','qb','mw','spectrum',cID{:},'Rewrite','yes');

end
