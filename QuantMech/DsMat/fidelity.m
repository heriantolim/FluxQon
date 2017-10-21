function F=fidelity(p1,p2)
%% Fidelity
%  F=fidelity(p1,p2) computes the fidelity of density matrix p1 and p2. If p1
%  and p2 are not of the same size, then the smaller matrix will be expanded to
%  the larger one by concatenation with zeros to the bottom and to the right.
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 03/08/2017
% Last modified: 03/08/2017

assert(iscomplexmatrix(p1) && iscomplexmatrix(p2),...
	'QuantMech:DsMat:fidelity:InvalidInput',...
	'Input to the density matrices must be a complex matrix.');
d1=size(p1,1);
d2=size(p2,1);

if d1<d2
	p1(d2,d2)=0;
elseif d2<d1
	p2(d1,d1)=0;
end

if isdiag(p1)
	p1=sqrt(p1);
else
	p1=sqrtm(p1);
end
p2=p1*p2*p1;
if isdiag(p2)
	p2=sqrt(p2);
else
	p2=sqrtm(p2);
end
F=abs(trace(p2));

end