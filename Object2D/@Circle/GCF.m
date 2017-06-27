function G=GCF(obj,x,y,z)
%% Geometrical Coupling Factor
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 24/06/2017
% Last modified: 24/06/2017

R=obj.Radius;
r=x.^2+y.^2;
Gr=R^2+r+z.^2;
Gz=R^2-r-z.^2;
r=sqrt(r);
a=Gr-2*R*r;
b=Gr+2*R*r;
[K,E]=ellipke(1-a./b);
b=sqrt(b);
Gz=2./a./b.*(Gz*E+a*K);
if r==0
	G=[0,0,Gz];
else
	Gr=2./a./b.*(Gr*E-a*K).*z./r;
	G=[Gr.*x./r,Gr.*y./r,Gz];
end

end