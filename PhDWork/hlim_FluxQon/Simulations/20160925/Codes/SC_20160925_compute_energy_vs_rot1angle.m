function SC_20160925_compute_energy_vs_rot1angle(a,varargin)

P=inputParser;
P.addRequired('Rot1Angle',@isrealvector);
P.addOptional('Rot2Angle',90,@isrealscalar);
P.addOptional('Rot3Angle',0,@isrealscalar);
P.addOptional('MagneticField',.08,@isrealscalar);
P.parse(a,varargin{:});
P=P.Results;

S=struct();
S.Rot1Angle=a;
S.Energy=cell(2,2,2);
a=a*pi/180;
b=P.Rot2Angle*pi/180;
c=P.Rot3Angle*pi/180;
N=numel(a);

ion=YSOEr();
ion.MagneticField=P.MagneticField;
for k=1:2
	for j=1:2
		for i=1:2
			ion.Isotope=i;
			ion.Site=j;
			ion.Class=k;
			S.Energy{i,j,k}=zeros(ion.HilbertDimension,N);
			for n=1:N
				ion.Rotation=[a(n),b,c];
				S.Energy{i,j,k}(:,n)=diag(ion.Hamiltonian);
			end
		end
	end
end
savedata(S,'R','energy','vs','rot1angle',...
	'in',P.MagneticField,'at',P.Rot2Angle,P.Rot3Angle,'Rewrite','yes');

end