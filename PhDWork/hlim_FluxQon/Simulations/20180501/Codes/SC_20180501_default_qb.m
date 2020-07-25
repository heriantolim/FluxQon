function qb=SC_20180501_default_qb(varargin)

P=inputParser;
P.addParameter('Frequency',2*pi*1.95e9,@isrealscalar);
P.addParameter('TunnelingFrequency',6.375814463646857e+12,@isrealscalar);
P.addParameter('Decoherence','on',@(x)any(strcmpi(x,{'on','off'})));
P.addParameter('Radius',5e-6,@isrealscalar);
P.parse(varargin{:});
P=P.Results;

qb=Circle3JJ(P.Frequency,P.TunnelingFrequency);
qb.Radius=P.Radius;

if strcmpi(P.Decoherence,'on')
	qb.DecayRate=2*pi*1e5;
	qb.DephasingRate=2*qb.DecayRate;
	qb.Temperature=2e-2;
end

end