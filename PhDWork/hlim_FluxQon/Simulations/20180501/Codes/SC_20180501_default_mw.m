function mw=SC_20180501_default_mw(varargin)

P=inputParser;
P.addParameter('Frequency',2*pi*1.95e9,@isrealscalar);
P.addParameter('FockDimension',3,@isintegerscalar);
P.addParameter('Decoherence','on',@(x)any(strcmpi(x,{'on','off'})));
P.addParameter('LossRate',2*pi*1e4,@isrealscalar);
P.parse(varargin{:});
P=P.Results;

mw=Microwave(P.Frequency,P.FockDimension);

if strcmpi(P.Decoherence,'on')
	mw.DecayRate=P.LossRate;
	mw.Temperature=2e-2;
end

end