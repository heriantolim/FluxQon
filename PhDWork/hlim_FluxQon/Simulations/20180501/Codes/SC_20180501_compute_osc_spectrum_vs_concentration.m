function cID=SC_20180501_compute_osc_spectrum_vs_concentration(n,varargin)

P=inputParser;
P.addRequired('NumPhotons',@isintegervector);
P.addOptional('Concentration',.005,@isrealvector);% atomic percent
P.addOptional('CalcMode',0,@isintegerscalar);
P.addOptional('Detuning',.5,@isrealscalar);% times Frequency
P.addOptional('CouplingRatio',1,@isrealscalar);
P.addOptional('ComponentID',{},@(x)isemptycell(x) || FileDir.iscomponentid(x));
P.addParameter('Frequency',2*pi*1.95e9,@isrealscalar);
P.addParameter('MaxProbeFreq',.1,@isrealscalar);% times Frequency
P.addParameter('TimeLength',1e-3,@isrealscalar);
P.parse(n,varargin{:});
P=P.Results;
cID=FileDir.expandlist(P.ComponentID);
cID=[{P.CalcMode,P.Detuning,P.CouplingRatio},cID];

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
	'Frequency',P.Frequency,'Decoherence',P.Decoherence);
Grms=runfunction({'S','C',20170621,'get','Grms'},ion,qb);
ion.CouplingStrength=Grms;

mw=runfunction({'default','mw'},...
	'Frequency',P.Frequency,'Decoherence',P.Decoherence);
MA=P.CouplingRatio*Grms*Constant.FluxQuantum ...
	/(pi*qb.Area*qb.TunnelingFrequency);

P.MaxProbeFreq=P.MaxProbeFreq*P.Frequency;
tSpan=0:pi/P.MaxProbeFreq:P.TimeLength;

S=struct();
S.Concentration=P.Concentration;
N=2^nextpow2(numel(tSpan));
S.Frequency=(0:2/N:1)*P.MaxProbeFreq;
i=0;
for n=P.NumPhotons
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
	psi0=State.kron(d,psi0{:});

	% The observables.
	O=Operator.kron(d,3,mw.Number);
	
	S.Intensity=cell(numel(S.Concentration),1);
	j=0;
	for CC=S.Concentration
		j=j+1;
		
		ion.NumIons=round(CC*ion.Volume/200 ...
			*Constant.AvogadroNumber*ion.Density/ion.MolarMass);
		
		mw.MagneticAmplitude=sqrt(ion.NumIons)*MA;
		
		if P.CalcMode<2
			V=Solve.UTE(tSpan,psi0,...
				Construct.Hamiltonian(qb,ion,mw,P.Interaction{:}),'Observable',O);
		else
			V=Solve.LME(tSpan,psi0*psi0',...
				Construct.Lindblad(qb,ion,mw,P.Interaction{:}),'Observable',O);	
		end
		
		% Fourier transform.
		V=fft(V,N);
		V=abs(V(1:N/2+1))/N;
		V(2:end-1)=2*V(2:end-1);
		
		S.Intensity{j}=V;
	end
	S.Intensity=cell2mat(S.Intensity);
	
	savedata(S,'R','osc','spectrum','vs','concentration',n,cID{:},'Rewrite','yes');
end

end