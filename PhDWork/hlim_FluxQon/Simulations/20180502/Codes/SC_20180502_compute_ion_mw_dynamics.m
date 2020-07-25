function cID=SC_20180502_compute_ion_mw_dynamics(varargin)

P=inputParser;
P.addRequired('SaveState',@isintegerscalar);
P.addRequired('NumPhotons',@isintegerscalar);
P.addOptional('CalcMode',0,@isintegerscalar);
P.addOptional('Concentration',.005,@isrealscalar);% atomic percent
P.addOptional('ComponentID',{},@(x)isemptycell(x) || FileDir.iscomponentid(x));
P.addParameter('Frequency',2*pi*1.95e9,@isrealscalar);
P.addParameter('TimeSpan',0:1e-10:5e-7,@isrealvector);
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

ion=cell(1,2);
for i=1:2
	ion{i}=runfunction({'S','C',20180501,'default','ion'},i,...
		'Frequency',P.Frequency,'FockDimension',n+2,...
		'Decoherence',P.Decoherence,'Concentration',P.Concentration);
end

mw=runfunction({'S','C',20180501,'default','mw'},...
	'Frequency',P.Frequency,'FockDimension',n+2,'Decoherence',P.Decoherence);

d=Hilbert.dimension(ion{:},mw);
M=numel(d);

% The initial state.
psi0=cell(2,M);
for i=1:M
	psi0{1,i}=i;% state index
	psi0{2,i}=eye(d(i),1);% ground state
end
i=3;
psi0{2,i}([1,n+1])=psi0{2,i}([n+1,1]);% excite the microwave state
psi0=State.kron(d,psi0{:});

% The observables.
O=cell(1,M);
O{1}=Operator.kron(d,1,ion{1}.Number(2));
O{2}=Operator.kron(d,2,ion{2}.Number(2));
O{3}=Operator.kron(d,3,mw.Number);

% Solve the evolution.
S=struct();
if P.SaveState==0
	if P.CalcMode<2
		[S.Time,S.Number]=Solve.UTE(P.TimeSpan,psi0,...
			Construct.Hamiltonian(ion{:},mw,P.Interaction{:}),...
			'Observable',O,'DispProgress',P.DispProgress);
	else
		[S.Time,S.Number]=Solve.LME(P.TimeSpan,psi0*psi0',...
			Construct.Lindblad(ion{:},mw,P.Interaction{:}),...
			'Observable',O,'DispProgress',P.DispProgress);
	end
else
	S.Dimension=d;
	if P.CalcMode<2
		[S.Time,S.State,S.Number]=Solve.UTE(P.TimeSpan,psi0,...
			Construct.Hamiltonian(ion{:},mw,P.Interaction{:}),...
			'Observable',O,'DispProgress',P.DispProgress);
	else
		[S.Time,S.State,S.Number]=Solve.LME(P.TimeSpan,psi0*psi0',...
			Construct.Lindblad(ion{:},mw,P.Interaction{:}),...
			'Observable',O,'DispProgress',P.DispProgress);
	end
end

% Saving ...
cID=[{'ion','mw','dynamics',n,P.CalcMode,P.Concentration},cID];
savedata(S,'R',cID,'Rewrite','yes');

end