function cID=SC_20180501_compute_fidelity_vs_g_ratio(n,varargin)

P=inputParser;
P.addRequired('NumPhotons',@isintegervector);
P.addOptional('CouplingRatio',1,@isrealvector);
P.addOptional('CalcMode',0,@isintegerscalar);
P.addOptional('Detuning',0,@isrealscalar);% times Frequency
P.addOptional('Concentration',.005,@isrealscalar);% atomic percent
P.addOptional('ComponentID',{},@(x)isemptycell(x) || FileDir.iscomponentid(x));
P.addParameter('Frequency',2*pi*1.95e9,@isrealscalar);
P.addParameter('TimeSpan',0:1e-11:6e-7);
P.parse(n,varargin{:});
P=P.Results;
cID=FileDir.expandlist(P.ComponentID);
tSpan=P.TimeSpan;

if isa(tSpan,'function_handle')
	tf=true;
elseif isrealvector(tSpan)
	tf=false;
else
	error('Invalid input to the parameter TimeSpan.');
end

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

qb=runfunction({'default','qb'},...
	'Frequency',P.Frequency*(1+P.Detuning),'Decoherence',P.Decoherence);

ion=runfunction({'default','ion'},1,...
	'Frequency',P.Frequency,'Decoherence',P.Decoherence,...
	'Concentration',P.Concentration);
Grms=runfunction({'S','C',20170621,'get','Grms'},ion,qb);
ion.CouplingStrength=Grms;

mw=runfunction({'default','mw'},...
	'Frequency',P.Frequency,'Decoherence',P.Decoherence);
MA=Grms*sqrt(ion.NumIons)*Constant.FluxQuantum ...
	/(pi*qb.Area*qb.TunnelingFrequency);

S=struct();
S.NumPhotons=P.NumPhotons;
S.CouplingRatio=P.CouplingRatio;
V=cell(numel(S.NumPhotons),numel(S.CouplingRatio));
S.Time=V;
S.Number=V;
S.Fidelity=V;
i=0;
for n=S.NumPhotons
	i=i+1;
	
	ion.FockDimension=n+2;
	mw.FockDimension=n+2;
	
	d=Hilbert.dimension(qb,ion,mw);
	M=numel(d);

	% The initial state.
	psi0=cell(2,M);
	for k=1:M
		psi0{1,k}=k;% state index
		psi0{2,k}=eye(d(k),1);% ground state
	end
	k=3;
	psi0{2,k}([1,n+1])=psi0{2,k}([n+1,1]);% excite the microwave state
	psi0_mw=psi0{2,k};% initial state of the microwave
	psi0=State.kron(d,psi0{:});

	% The observables.
	O=cell(1,2);
	O{1}=Operator.kron(d,1,qb.Number);
	O{2}=@(x)fidelity(psi0_mw,subdsmat(x,d,2));
	
	j=0;
	for CR=S.CouplingRatio
		j=j+1;
		if tf
			tSpan=P.TimeSpan(CR);
		end
		
		mw.MagneticAmplitude=CR*MA;
		
		if P.CalcMode<2
			[S.Time{i,j},V]=Solve.UTE(tSpan,psi0,...
				Construct.Hamiltonian(qb,ion,mw,P.Interaction{:}),'Observable',O);
		else
			[S.Time{i,j},V]=Solve.LME(tSpan,psi0*psi0',...
				Construct.Lindblad(qb,ion,mw,P.Interaction{:}),'Observable',O);	
		end
		
		S.Number(i,j)=V(1);
		S.Fidelity(i,j)=V(2);
	end
end

% Saving ...
cID=[{'fidelity','vs','g','ratio',...
	P.CalcMode,P.Detuning,P.Concentration},cID];
savedata(S,'R',cID,'Rewrite','yes');

end