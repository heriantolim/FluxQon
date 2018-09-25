function Gr=GCFr(R,r,z)

if r==0
	Gr=0;
else
	Gr=R^2+r.^2+z.^2;
	a=Gr-2*R*r;
	b=Gr+2*R*r;
	[K,E]=ellipke(1-a./b);
	Gr=2./a./sqrt(b).*(Gr.*E-a.*K).*z./r;
end

end