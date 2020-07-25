function SC_20160715_compute_frequency_vs_tuning(PT,varargin)

P=inputParser;
P.addRequired('TuningPhase',@isrealvector);
P.addOptional('JosephsonEnergy',5.203156259013710e-22,@isrealscalar);
P.addOptional('AlphaCoefficient',.7,@isrealscalar);
P.parse(PT,varargin{:});
P=P.Results;

M=100;
qb=a3JJ(P.JosephsonEnergy/10,P.JosephsonEnergy,P.AlphaCoefficient);

N=numel(PT);
w=zeros(2,N);

fprintf('Progress: ');
for j=1:N
	fprintf('%2d%%',floor(100*j/N));
	qb.PhaseTune=PT(j);
	qb.solveTISE(M);
	w(1,j)=qb.Frequency;
	w(2,j)=abs(qb.TunnelingFrequency);
	fprintf('\b\b\b');
end
fprintf(repmat('\b',1,11));
savedata([PT;w/min(w(1,:))],'R','frequency','vs','tuning','Rewrite','yes');

end