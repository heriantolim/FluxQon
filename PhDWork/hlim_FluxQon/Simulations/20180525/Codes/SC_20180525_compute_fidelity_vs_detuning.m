function cID=SC_20180525_compute_fidelity_vs_detuning(varargin)

P=inputParser;
P.addOptional('Detuning',0,@isrealvector);% times g_eff
P.addOptional('CalcMode',0,@isintegerscalar);
P.addOptional('NumPhotons',1,@isintegerscalar);
P.addOptional('NumIons',1,@isrealscalar);
P.addOptional('CouplingRatio',1,@isrealscalar);
P.addOptional('LineStrength3',1e3,@isrealscalar);
P.addOptional('LossRate',[1e6,1e8],...% [mw,op]
	@(x)isrealvector(x) && numel(x)==2);
P.addOptional('ComponentID',{},@(x)isemptycell(x) || FileDir.iscomponentid(x));
P.addParameter('MwFrequency',2*pi*1.95e9,@isrealscalar);
P.addParameter('TimeSpan',0:1e-4:1,@isrealvector);
P.parse(varargin{:});
P=P.Results;
n=P.NumPhotons;
N=P.NumIons;
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

g_eff=2*pi*P.LineStrength3*P.CouplingRatio;

qb=runfunction({'S','C',20180501,'default','qb'},'Decoherence',P.Decoherence);

ion=cell(1,N);
for j=1:N
	ion{j}=runfunction({'default','ion'},1,...
		'MwFrequency',P.MwFrequency,'Decoherence',P.Decoherence,...
		'LineStrength',2*pi*P.LineStrength3*[0,10,1]);
end
G=runfunction({'S','C',20170621,'get','Grms'},ion{1},qb);
G(2,1)=g_eff;
for j=1:N
	ion{j}.CouplingStrength=G;
end

mw=runfunction({'S','C',20180501,'default','mw'},...
	'Frequency',P.MwFrequency,'FockDimension',n+2,...
	'Decoherence',P.Decoherence,'LossRate',2*pi*P.LossRate(1));
mw.MagneticAmplitude=g_eff ...
	*Constant.FluxQuantum/(pi*qb.Area*qb.TunnelingFrequency);

op=runfunction({'S','C',20180521,'default','op'},...
	'Frequency',diff(ion{1}.Frequency([1,4])),'FockDimension',n+2,...
	'Decoherence',P.Decoherence,'LossRate',2*pi*P.LossRate(2));

d=Hilbert.dimension(qb,ion{:},mw,op);
M=numel(d);

% The initial state.
psi0=cell(2,M);
for i=1:M
	psi0{1,i}=i;% state index
	psi0{2,i}=eye(d(i),1);% ground state
end
i=N+3;
psi0{2,i}([1,n+1])=psi0{2,i}([n+1,1]);% excite the optical state
psi0_op=psi0{2,i};% initial state of the optical
psi0=State.kron(d,psi0{:});

% The observables.
O=@(x)fidelity(psi0_op,subdsmat(x,d,N+2));

S=struct();
S.Detuning=P.Detuning;
J=numel(S.Detuning);
S.Time=cell(1,J);
S.Fidelity=cell(1,J);

for j=1:J
	tSpan=P.TimeSpan*sqrt(n/N/P.CouplingRatio);
	qb.Frequency=P.MwFrequency+S.Detuning(j)*g_eff;
	
	% Solve the evolution.
	if P.CalcMode<2
		[S.Time{j},S.Fidelity{j}]=Solve.UTE(tSpan,psi0,...
			Construct.Hamiltonian(qb,ion{:},mw,op,P.Interaction{:}),'Observable',O);
	else
		[S.Time{j},S.Fidelity{j}]=Solve.LME(tSpan,psi0*psi0',...
			Construct.Lindblad(qb,ion{:},mw,op,P.Interaction{:}),'Observable',O);	
	end
end

% Saving ...
if P.CalcMode>=2
	cID=[num2cell(P.LossRate),cID];
end
cID=[{'fidelity','vs','detuning',...
	P.CalcMode,P.NumPhotons,P.NumIons,P.CouplingRatio,P.LineStrength3},cID];
savedata(S,'R',cID,'Rewrite','yes');

end