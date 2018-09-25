function Gz=GCFz(R,r,z)

Gr=R^2+r.^2+z.^2;
Gz=R^2-r.^2-z.^2;
a=Gr-2*R*r;
b=Gr+2*R*r;
[K,E]=ellipke(1-a./b);
Gz=2./a./sqrt(b).*(Gz.*E+a.*K);

end