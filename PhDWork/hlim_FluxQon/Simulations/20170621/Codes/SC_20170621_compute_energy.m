function SC_20170621_compute_energy(A,b,c)

P=inputParser;
P.addRequired('Rot1Angle',@isrealvector);
P.addRequired('Rot2Angle',@isrealscalar);
P.addRequired('Rot3Angle',@isrealscalar);
P.parse(A,b,c);

Na=numel(A);
dE=[reshape(A,1,Na);zeros(4,Na)];
ion=YSOEr(...
	'Multiplet',[1,2],...
	'MagneticField',1,...
	'RoundRelTol',0);
for s=1:2
	ion.Site=s;
	for h=1:Na
		ion.Rotation=[A(h),b,c]*pi/180;
		E=diag(ion.Hamiltonian);
		dE(1+s,h)=E(2)-E(1);
		dE(3+s,h)=E(4)-E(3);
	end
end
savedata(dE,'R','energy',b,c,'Rewrite','yes');

end