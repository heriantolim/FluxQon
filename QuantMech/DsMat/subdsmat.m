function pn=subdsmat(p,d,n)
%% Sub-system Density Matrix 
%  pn=subdsmat(p,d,n) computes the density matrix of the subsystem of index n
%  from the density matrix p, by performing partial trace of p over all other
%  sub-systems. This function is a complementary operation to the partial trace
%  (ptrace). The integer vector d is the list of the subspace dimensions.
%
% See also: ptrace.
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 03/08/2017
% Last modified: 03/08/2017

assert(iscomplexmatrix(p),...
	'QuantMech:DsMat:subdsmat:InvalidInput',...
	'Input to the density matrix must be a complex matrix.');
assert(isintegervector(d),...
	'QuantMech:DsMat:subdsmat:InvalidInput',...
	'Input to the subspace dimensions must be a positive integer vector.');
D=prod(d);
assert(D==size(p,1),...
	'QuantMech:DsMat:subdsmat:InvalidInput',...
	['The product of the elements in the subspace dimension must be ',...
		'equal to the dimension of the density matrix.']);
N=numel(d);
assert(isintegerscalar(n) && n>0 && n<=N,...
	'QuantMech:DsMat:subdsmat:InvalidInput',...
	['Input to the subspace index must be an integer scalar between 1 ',...
		'and the total number of the subspaces.']);

if N==1
	pn=p;
	return
end
d=[prod(d(1:(n-1))),d(n),prod(d((n+1):N))];

if d(1)>0
	D=d(2)*d(3);
	pn=zeros(D);
	for k=1:D
		for l=1:D
			for i=1:d(1)
				pn(k,l)=pn(k,l)+p((i-1)*D+k,(i-1)*D+l);
			end
		end
	end
	p=pn;
end

if d(3)>0
	D=d(2);
	pn=zeros(D);
	for i=1:D
		for j=1:D
			for k=1:d(3)
				pn(i,j)=pn(i,j)+p((i-1)*d(3)+k,(j-1)*d(3)+k);
			end
		end
	end
end

end