function E=SC_20170621_get_frequency(m,varargin)

P=inputParser;
if isnumeric(m)
	P.addRequired('Multiplet',@isintegervector);
	P.addRequired('Site',@isintegervector);
	P.addRequired('Rot1Angle',@isrealvector);
	P.addOptional('Rot2Angle',270,@isrealscalar);
	P.addOptional('Rot3Angle',0,@isrealscalar);
	P.addOptional('MagneticField',1,@isrealscalar);
	P.parse(m,varargin{:});
	P=P.Results;
	
	m=P.Multiplet;
	s=P.Site;
	a=P.Rot1Angle;
	b=P.Rot2Angle;
	c=P.Rot3Angle;
	B=P.MagneticField;
else
	P.addRequired('Ion',@(x)isa(x,'Ion') && isa(x,'Cylinder'));
	P.parse(m,varargin{:});
	P=P.Results;
	
	m=P.Ion.Multiplet;
	s=P.Ion.Site;
	a=P.Ion.Rotation*180/pi;
	b=a(2); c=a(3); a=a(1);
	B=P.Ion.MagneticField(3);
end

E=loaddata('R','energy',b,c);
I=2*mod(m(:)-1,2)+mod(s(:)'-1,2)+2;
J=arrayfun(@(x)find(x==E(1,:)),mod(a,180));

if isvector(I)
	E=reshape(E(I,J),numel(I),numel(J));
else
	E=reshape(E(I(:),J),[size(I),numel(J)]);
end

E=E*B/Constant.ReducedPlanck;

end