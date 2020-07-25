function ion=SC_20180501_default_ion(varargin)

P=inputParser;
P.addOptional('Site',1,@isintegerscalar);
P.addParameter('Frequency',2*pi*1.95e9,@isrealscalar);
P.addParameter('FockDimension',3,@isintegerscalar);
P.addParameter('Decoherence','on',@(x)any(strcmpi(x,{'on','off'})));
P.addParameter('Concentration',.005,@isrealscalar);% atomic percent
P.addParameter('LineStrength',2*pi*10,@isrealscalar);
P.addParameter('Radius',4e-6,@isrealscalar);
P.addParameter('Height',8e-6,@isrealscalar);
P.parse(varargin{:});
P=P.Results;

a=152;
B=runfunction({'S','C',20170621,'get','frequency'},1,1,a,270,0);
B=P.Frequency/B;

ion=CylinderYSOEr(...
	'FockDimension',P.FockDimension,...
	'Isotope',2,...
	'Multiplet',1,...
	'Site',P.Site,...
	'Class',1,...
	'Rotation',[a,270,0]*pi/180,...
	'MagneticField',B,...
	'Radius',P.Radius,...
	'Height',P.Height,...
	'LineStrength',P.LineStrength);

% NumIons = 1/2 * Concentration * Volume * NumberDensity
% half because there are two types of ions in the crystal.
ion.NumIons=round(P.Concentration*ion.Volume/200 ...
	*Constant.AvogadroNumber*ion.Density/ion.MolarMass);

if strcmpi(P.Decoherence,'on')
	ion.DecayRate=2*pi*1e6;
	ion.Temperature=2e-2;
end

end
