function pn=ptrace(p,d,n)
%% Partial Trace
%  pn=ptrace(p,d,n) computes the partial trace of the density matrix, p, over
%  the subsystem(s) of index n. The input to p can also be a state (a complex
%  vector), in which case p will be parsed as p*p'. The integer vector d is the
%  list of the subspace dimensions. The index n can be inputted as an integer
%  scalar or vector. The elements in n must not have duplicates. This function
%  is a complementary operation to subdsmat.
%
% See also: subdsmat.
%
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 03/08/2017
% Last modified: 27/04/2018

D=size(p);
tf=any(D==1);
assert(iscomplexmatrix(p) && (tf || D(1)==D(2)),...
	'FluxQon:State:ptrace:InvalidInput',...
	'The first input must be a vector or a square matrix of complex numbers.');
if tf
	if D(1)==1
		p=p.'*conj(p);
	else
		p=p*p';
	end
end
D=max(D);

assert(isintegervector(d),...
	'FluxQon:State:ptrace:InvalidInput',...
	'Input to the subspace dimensions must be a positive integer vector.');
assert(D==prod(d),...
	'FluxQon:State:ptrace:InvalidInput',...
	['The product of the elements in the subspace dimension must be ',...
		'equal to the dimension of the density matrix.']);
N=numel(d);

assert(isintegervector(n) && all(n>0 & n<=N) && numel(unique(n))==numel(n),...
	'FluxQon:State:ptrace:InvalidInput',...
	['Input to the subspace index must be an integer scalar or vector, ',...
		'with no duplicates, between 1 and the total number of the subspaces.']);

if N==1
	pn=trace(p);
	return
end

while true
	q=[prod(d(1:(n(1)-1))),d(n(1)),prod(d((n(1)+1):N))];
	r=q(2)*q(3);
	pn=zeros(q(1)*q(3));
	for i=1:q(1)
		for j=1:q(1)
			for m=1:q(2)
				for k=1:q(3)
					for l=1:q(3)
						pn((i-1)*q(3)+k,(j-1)*q(3)+l)=pn((i-1)*q(3)+k,(j-1)*q(3)+l) ...
							+p((i-1)*r+(m-1)*q(3)+k,(j-1)*r+(m-1)*q(3)+l);
					end
				end
			end
		end
	end
	
	if numel(n)>1
		p=pn;
		d(n(1))=[];
		N=N-1;
		ix=n>n(1);
		n(ix)=n(ix)-1;
		n(1)=[];
	else
		break
	end
end

end