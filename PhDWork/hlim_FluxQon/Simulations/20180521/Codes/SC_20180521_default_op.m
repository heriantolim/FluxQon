function op=SC_20180521_default_op(varargin)

P=inputParser;
P.addParameter('Frequency',2*pi*1.95e14,@isrealscalar);
P.addParameter('FockDimension',3,@isintegerscalar);
P.addParameter('Decoherence','on',@(x)any(strcmpi(x,{'on','off'})));
P.addParameter('LossRate',2*pi*1e8,@isrealscalar);
P.parse(varargin{:});
P=P.Results;

op=Optical(P.Frequency,P.FockDimension);

if strcmpi(P.Decoherence,'on')
	op.DecayRate=P.LossRate;
	op.Temperature=2e-2;
end

end