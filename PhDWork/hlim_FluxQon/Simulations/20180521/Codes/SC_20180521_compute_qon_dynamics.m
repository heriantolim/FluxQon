function cID=SC_20180521_compute_qon_dynamics(varargin)

P=inputParser;
P.addRequired('SaveState',@isintegerscalar);
P.addOptional('CalcMode',0,@isintegerscalar);
P.addOptional('NumPhotons',1,@isintegerscalar);
P.addOptional('NumIons',1e4,@isrealscalar);
P.addOptional('Site',1,@isintegerscalar);
P.addOptional('QbDetuning',0,@isrealscalar);% times Frequency
P.addOptional('CouplingRatio',1e2,@isrealscalar);
P.addOptional('LineStrength3',1e3,@isrealscalar);
P.addOptional('LossRate',[1e6,1e8],...% [mw,op]
	@(x)isrealvector(x) && numel(x)==2);
P.addOptional('ComponentID',{},@(x)isemptycell(x) || FileDir.iscomponentid(x));
P.addParameter('MwFrequency',2*pi*1.95e9,@isrealscalar);
P.addParameter('TimeSpan',0:1e-7:1e-4,@isrealvector);
P.addParameter('DispProgress',true);
P.parse(varargin{:});
P=P.Results;
n=P.NumPhotons;
cID=FileDir.expandlist(P.ComponentID);

switch P.CalcMode
	case {0,1}
		P.Decoherence='off';
	case {2,3}
		P.Decoherence='on';
end

switch P.CalcMode
	case {0,2}
		P.Interaction={};
	case {1,3}
		P.Interaction={'rwa'};
end

qb=runfunction({'S','C',20180501,'default','qb'},...
	'Frequency',P.MwFrequency*(1+P.QbDetuning),...
	'Decoherence',P.Decoherence);

ion=runfunction({'default','ion'},1);
P.OpFrequency=diff(ion.Frequency([1,4]));
G=round(ion.Rotation/pi*180);
G=runfunction({'S','C',20170621,'get','Grms'},1,1:2,G(1),G(2),G(3),...
	ion.Radius,ion.Height,qb.Radius,qb.TunnelingFrequency).';
G=[G;2*pi*P.LineStrength3*P.CouplingRatio*[1,1]];

ion=runfunction({'default','ion'},P.NumIons,P.Site,...
	'MwFrequency',P.MwFrequency,'FockDimension',n+2,...
	'Decoherence',P.Decoherence,...
	'CouplingStrength',G(:,P.Site),...
	'LineStrength',2*pi*P.LineStrength3*[0,10,1]);

mw=runfunction({'S','C',20180501,'default','mw'},...
	'Frequency',P.MwFrequency,'FockDimension',n+2,...
	'Decoherence',P.Decoherence,'LossRate',2*pi*P.LossRate(1));
mw.MagneticAmplitude=G(2,1) ...
	*Constant.FluxQuantum/(pi*qb.Area*qb.TunnelingFrequency);

op=runfunction({'default','op'},...
	'Frequency',P.OpFrequency,'FockDimension',n+2,...
	'Decoherence',P.Decoherence,'LossRate',2*pi*P.LossRate(2));

d=Hilbert.dimension(qb,ion,mw,op);
M=numel(d);

% The initial state.
psi0=cell(2,M);
for i=1:M
	psi0{1,i}=i;% state index
	psi0{2,i}=eye(d(i),1);% ground state
end
i=4;
psi0{2,i}([1,n+1])=psi0{2,i}([n+1,1]);% excite the optical state
psi0_op=psi0{2,i};% initial state of the optical
psi0=State.kron(d,psi0{:});

% The observables.
O=cell(1,7);
O{1}=Operator.kron(d,1,qb.Number);
O{2}=Operator.kron(d,2,ion.Number(2));
O{3}=Operator.kron(d,2,ion.Number(3));
O{4}=Operator.kron(d,2,ion.Number(4));
O{5}=Operator.kron(d,3,mw.Number);
O{6}=Operator.kron(d,4,op.Number);
O{7}=@(x)fidelity(psi0_op,subdsmat(x,d,3));

% Solve the evolution.
S=struct();
if P.SaveState==0
	if P.CalcMode<2
		[S.Time,S.Number]=Solve.UTE(P.TimeSpan,psi0,...
			Construct.Hamiltonian(qb,ion,mw,op,P.Interaction{:}),...
			'Observable',O,'DispProgress',P.DispProgress);
	else
		[S.Time,S.Number]=Solve.LME(P.TimeSpan,psi0*psi0',...
			Construct.Lindblad(qb,ion,mw,op,P.Interaction{:}),...
			'Observable',O,'DispProgress',P.DispProgress);
	end
else
	S.Dimension=d;
	if P.CalcMode<2
		[S.Time,S.State,S.Number]=Solve.UTE(P.TimeSpan,psi0,...
			Construct.Hamiltonian(qb,ion,mw,op,P.Interaction{:}),...
			'Observable',O,'DispProgress',P.DispProgress);
	else
		[S.Time,S.State,S.Number]=Solve.LME(P.TimeSpan,psi0*psi0',...
			Construct.Lindblad(qb,ion,mw,op,P.Interaction{:}),...
			'Observable',O,'DispProgress',P.DispProgress);
	end
end
S.Fidelity=S.Number{end};
S.Number(end)=[];

% Saving ...
if P.CalcMode>=2
	cID=[num2cell(P.LossRate),cID];
end
cID=[{'qon','dynamics',P.CalcMode,n,P.NumIons,P.Site,...
	P.QbDetuning,P.CouplingRatio,P.LineStrength3},cID];
savedata(S,'R',cID,'Rewrite','yes');

end
