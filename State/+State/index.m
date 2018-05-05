function p=index(d,n,i)
%% State Index
%  p=State.index(d,n,i) returns the indices of a list of states i, belonging
%  to the n-th object, in a kronecker product of Hilbert spaces whose dimensions
%  are listed in the vector d.
%
% Inputs:
%  d: Subspace dimensions, specified as a positive integer vector.
%  n: Subspace indices, specified as a positive integer scalar.
%  i: List of the state indices being queried, specified as a positive integer
%     vector.
%
% See also: kron.
%
% Requires package:
%  - MatCommon_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
% 
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 19/12/2015
% Last modified: 19/12/2015

assert(isintegervector(d) && all(d>0),...
	'FluxQon:State:index:InvalidInput',...
	'Input to the subspace dimensions must be a positive integer vector.');
N=numel(d);
assert(isintegerscalar(n) && n>0 && n<=N,...
	'FluxQon:State:index:InvalidInput',...
	['Input to the subspace index must be an integer scalar between 1 ',...
		'and the total number of the subspaces.']);
assert(isintegervector(i) && all(i>0) && all(i<=d(n)),...
	'FluxQon:State:index:InvalidInput',...
	'Input to the state indices must be a positive integer vector.');

D1=prod(d(1:(n-1)));
D2=prod(d((n+1):N));
K=numel(i);
p=zeros(1,D1*D2*K);
j=0;
m=0;
while j<D1
	k=0;
	while k<K
		k=k+1;
		l=0;
		while l<D2
			l=l+1;
			m=m+1;
			p(m)=(j*d(n)+i(k)-1)*D2+l;
		end
	end
	j=j+1;
end

end