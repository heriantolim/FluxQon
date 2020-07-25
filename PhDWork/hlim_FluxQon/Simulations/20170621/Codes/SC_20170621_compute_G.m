function SC_20170621_compute_G(A,b,c)

P=inputParser;
P.addRequired('Rot1Angle',@isrealvector);
P.addRequired('Rot2Angle',@isrealscalar);
P.addRequired('Rot3Angle',@isrealscalar);
P.parse(A,b,c);

C=Circle(1);
ion=YSOEr(2);
ion.Class=1;
ion.MagneticField=1;
ion.RoundRelTol=0;
S=struct();
S.X=-.9:.01:.9;
S.Y=S.X;
S.Z=0:.01:1.8;
S.G=cell(2);
Nx=numel(S.X);
Ny=numel(S.Y);
Nz=numel(S.Z);

for a=A
	ion.Rotation=[a,b,c]*pi/180;
	for m=1:2
		ion.Multiplet=m;
		for s=1:2
			ion.Site=s;
			ion.Hamiltonian('S');
			S.G{m,s}=zeros(Nx,Ny,Nz);
			for i=1:Nx
				for j=1:Ny
					for k=1:Nz
						HZ=ion.ZeemanHamiltonian(Rotation.Y(-pi/2) ...
							*C.GCF(S.X(i),S.Y(j),S.Z(k)));
						S.G{m,s}(i,j,k)=Constant.VacuumPermeability ...
							/4/Constant.FluxQuantum*abs(HZ(1,2));
					end
				end
			end
		end
	end
	savedata(S,'R','G',a,b,c,'Rewrite','yes');
end

end