function G=SC_20170621_get_Grms(m,varargin)

P=inputParser;
if isnumeric(m)
	P.addRequired('Multiplet',@isintegervector);
	P.addRequired('Site',@isintegervector);
	P.addRequired('Rot1Angle',@isrealvector);
	P.addOptional('Rot2Angle',270,@isrealscalar);
	P.addOptional('Rot3Angle',0,@isrealscalar);
	P.addOptional('CrystalRadius',.8,@isrealscalar);
	P.addOptional('CrystalHeight',1.6,@isrealscalar);
	P.addOptional('QubitRadius',1,@isrealscalar);
	P.addOptional('TunnelingFrequency',1,@isrealscalar);
	P.parse(m,varargin{:});
	P=P.Results;
	
	m=P.Multiplet;
	s=P.Site;
	a=P.Rot1Angle;
	b=P.Rot2Angle;
	c=P.Rot3Angle;
	R=P.QubitRadius;
	wx=P.TunnelingFrequency;
	r=P.CrystalRadius/R;
	z=P.CrystalHeight/R;
else
	P.addRequired('Ion',@(x)isa(x,'Ion') && isa(x,'Cylinder'));
	P.addOptional('Qubit',[],@(x)isa(x,'Qubit') && isa(x,'Circle'));
	P.parse(m,varargin{:});
	P=P.Results;
	
	m=P.Ion.Multiplet;
	s=P.Ion.Site;
	a=P.Ion.Rotation*180/pi;
	b=a(2); c=a(3); a=a(1);
	if isempty(P.Qubit)
		R=1;
		wx=1;
		r=.8;
		z=1.6;
	else
		R=P.Qubit.Radius;
		wx=P.Qubit.TunnelingFrequency;
		r=P.Ion.Radius/R;
		z=P.Ion.Height/R;
	end
end

G=loaddata('R','Grms',b,c,r,z);
I=2*mod(m(:)-1,2)+mod(s(:)'-1,2)+2;
J=arrayfun(@(x)find(x==G(1,:)),mod(a,180));

if isvector(I)
	G=reshape(G(I,J),numel(I),numel(J))/R*wx;
else
	G=reshape(G(I(:),J),[size(I),numel(J)])/R*wx;
end

end