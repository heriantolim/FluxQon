function SC_20160625_compute_energy_near_pi(PB,varargin)

P=inputParser;
P.addRequired('BiasPhase',@isrealvector);
P.addOptional('JosephsonEnergy',Constant.FluxQuantum/2/pi*1e-6,@isrealscalar);
P.parse(PB,varargin{:});
P=P.Results;

M=1000;
qb=SQUID(P.JosephsonEnergy/25,P.JosephsonEnergy*2/pi,P.JosephsonEnergy);

N=numel(PB);
EN=zeros(10,N);

fprintf('Progress: ');
for j=1:N
	fprintf('%2d%%',floor(100*j/N));
	qb.PhaseBias=PB(j);
	[E,~]=qb.solveTISE(M);
	EN(:,j)=E(1:10).';
	fprintf('\b\b\b');
end
fprintf(repmat('\b',1,11));
savedata([PB;EN/P.JosephsonEnergy],'R','energy','near','pi','Rewrite','yes');

end