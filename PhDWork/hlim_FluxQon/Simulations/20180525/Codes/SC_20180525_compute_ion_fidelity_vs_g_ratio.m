function cID=SC_20180525_compute_ion_fidelity_vs_g_ratio(varargin)

P=inputParser;
P.addOptional('CouplingRatio',1,@isrealvector);
P.addOptional('CalcMode',0,@isintegerscalar);
P.addOptional('NumPhotons',1,@isintegerscalar);
P.addOptional('NumIons',1,@isrealscalar);
P.addOptional('LineStrength3',1e3,@isrealscalar);
P.addOptional('LossRate',[1e6,1e8],...% [mw,op]
	@(x)isrealvector(x) && numel(x)==2);
P.addOptional('ComponentID',{},@(x)isemptycell(x) || FileDir.iscomponentid(x));
P.addParameter('MwFrequency',2*pi*1.95e9,@isrealscalar);
P.addParameter('TimeSpan',0:1e-5:1e-1,@isrealvector);
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

ion=runfunction({'default','ion'},1);
P.OpFrequency=diff(ion.Frequency([1,4]));

ion=cell(1,N);
for i=1:N
	ion{i}=runfunction({'default','ion'},1,...
		'MwFrequency',P.MwFrequency,'Decoherence',P.Decoherence);
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
O=@(x)fidelity(psi0_op,subdsmat(x,d,N+1));

S=struct();
S.CouplingRatio=P.CouplingRatio;
J=numel(S.CouplingRatio);
S.Time=cell(1,J);
S.Fidelity=cell(1,J);

ls=2*pi*P.LineStrength3*[10,1,1,10];
for j=1:J
	cr=S.CouplingRatio(j);
	tSpan=P.TimeSpan*sqrt(n/N/cr);
	for i=1:N
		ion{i}.LineStrength=[cr*ls(2),ls,ls(3)*cr];
	end
	
	% Solve the evolution.
	if P.CalcMode<2
		[S.Time{j},S.Fidelity{j}]=Solve.UTE(tSpan,psi0,...
			Construct.Hamiltonian(ion{:},mw,op,P.Interaction{:}),'Observable',O);
	else
		[S.Time{j},S.Fidelity{j}]=Solve.LME(tSpan,psi0*psi0',...
			Construct.Lindblad(ion{:},mw,op,P.Interaction{:}),'Observable',O);	
	end
end

% Saving ...
if P.CalcMode>=2
	cID=[num2cell(P.LossRate),cID];
end
cID=[{'ion','fidelity','vs','g','ratio',...
	P.CalcMode,P.NumPhotons,P.NumIons,P.LineStrength3},cID];
savedata(S,'R',cID,'Rewrite','yes');

end