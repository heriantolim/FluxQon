function SC_20170621_compute_Grms(A,b,c,r,z)

P=inputParser;
P.addRequired('Rot1Angle',@isrealvector);
P.addRequired('Rot2Angle',@isrealscalar);
P.addRequired('Rot3Angle',@isrealscalar);
P.addRequired('CrystalRadius',@isrealscalar);% as a fraction of the qubit radius
P.addRequired('CrystalHeight',@isrealscalar);% as a fraction of the qubit radius
P.parse(A,b,c,r,z);

r2=r^2;
Na=numel(A);
G=[reshape(A,1,Na);nan(4,Na)];
for h=1:Na
	g=1;
	S=loaddata('R','G',A(h),b,c);
	if ~isempty(S)
		Nx=numel(S.X);
		Ny=numel(S.Y);
		Nz=find(S.Z<=z,1,'last');
		for m=1:2
			for s=1:2
				g=g+1;
				G(g,h)=0;
				n=0;
				for i=1:Nx
					for j=1:Ny
						if S.X(i)^2+S.Y(j)^2<=r2
							n=n+Nz;
							G(g,h)=G(g,h)+sum(S.G{m,s}(i,j,1:Nz).^2);
						end
					end
				end
				G(g,h)=sqrt(G(g,h)/n);
			end
		end
	end
end
savedata(G,'R','Grms',b,c,r,z,'Rewrite','yes');

end