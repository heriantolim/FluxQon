function cID=SC_20180525_compute_ion_dynamics(varargin)

P=inputParser;
P.addRequired('SaveState',@isintegerscalar);
P.addOptional('CalcMode',0,@isintegerscalar);
P.addOptional('NumPhotons',1,@isintegerscalar);
P.addOptional('NumIons',1,@isintegerscalar);
P.addOptional('Site',1,@isintegerscalar);
P.addOptional('CouplingRatio',1,@isrealscalar);
P.addOptional('LineStrength3',1e3,@isrealscalar);
P.addOptional('LossRate',[1e6,1e8],...% [mw,op]
	@(x)isrealvector(x) && numel(x)==2);
P.addOptional('ComponentID',{},@(x)isemptycell(x) || FileDir.iscomponentid(x));
P.addParameter('MwFrequency',2*pi*1.95e9,@isrealscalar);
P.addParameter('TimeSpan',0:1e-5:1e-2,@isrealvector);
P.addParameter('DispProgress',true);
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

if P.Site==0
	s=mod((1:N)-1,2)+1;
else
	s=P.Site*ones(1,N);
end

ion=runfunction({'default','ion'},1);
P.OpFrequency=diff(ion.Frequency([1,4]));

ion=cell(1,N);
for j=1:N
	ion{j}=runfunction({'default','ion'},s(j),...
		'MwFrequency',P.MwFrequency,'Decoherence',P.Decoherence,...
		'LineStrength',2*pi*P.LineStrength3*[P.CouplingRatio,10,1]);
end

mw=runfunction({'S','C',20180501,'default','mw'},...
	'Frequency',P.MwFrequency,'FockDimension',n+2,...
	'Decoherence',P.Decoherence,'LossRate',2*pi*P.LossRate(1));

op=runfunction({'S','C',20180521,'default','op'},...
	'Frequency',P.OpFrequency,'FockDimension',n+2,...
	'Decoherence',P.Decoherence,'LossRate',2*pi*P.LossRate(2));

d=Hilbert.dimension(ion{:},mw,op);
M=numel(d);

% The initial state.
psi0=cell(2,M);
for i=1:M
	psi0{1,i}=i;% state index
	psi0{2,i}=eye(d(i),1);% ground state
end
i=N+2;
psi0{2,i}([1,n+1])=psi0{2,i}([n+1,1]);% excite the optical state
psi0_op=psi0{2,i};% initial state of the optical
psi0=State.kron(d,psi0{:});

% The observables.
O=cell(1,4*N+3);
i=0; k=0;
for j=1:N
	i=i+1;
	for l=1:4
		k=k+1; O{k}=Operator.kron(d,i,ion{j}.Number(l));
	end
end
i=i+1; k=k+1; O{k}=Operator.kron(d,i,mw.Number);
i=i+1; k=k+1; O{k}=Operator.kron(d,i,op.Number);
i=i-1; k=k+1; O{k}=@(x)fidelity(psi0_op,subdsmat(x,d,i));

% Solve the evolution.
S=struct();
if P.SaveState==0
	if P.CalcMode<2
		[S.Time,S.Number]=Solve.UTE(P.TimeSpan,psi0,...
			Construct.Hamiltonian(ion{:},mw,op,P.Interaction{:}),...
			'Observable',O,'DispProgress',P.DispProgress);
	else
		[S.Time,S.Number]=Solve.LME(P.TimeSpan,psi0*psi0',...
			Construct.Lindblad(ion{:},mw,op,P.Interaction{:}),...
			'Observable',O,'DispProgress',P.DispProgress);
	end
else
	S.Dimension=d;
	if P.CalcMode<2
		[S.Time,S.State,S.Number]=Solve.UTE(P.TimeSpan,psi0,...
			Construct.Hamiltonian(ion{:},mw,op,P.Interaction{:}),...
			'Observable',O,'DispProgress',P.DispProgress);
	else
		[S.Time,S.State,S.Number]=Solve.LME(P.TimeSpan,psi0*psi0',...
			Construct.Lindblad(ion{:},mw,op,P.Interaction{:}),...
			'Observable',O,'DispProgress',P.DispProgress);
	end
end
S.Fidelity=S.Number{end};
S.Number(end)=[];

% Saving ...
if P.CalcMode>=2
	cID=[num2cell(P.LossRate),cID];
end
cID=[{'ion','dynamics',P.CalcMode,n,P.NumIons,P.Site,...
	P.CouplingRatio,P.LineStrength3},cID];
savedata(S,'R',cID,'Rewrite','yes');

end
