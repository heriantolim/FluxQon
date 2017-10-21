function H=kron(d,varargin)
%% Kronecker Product of Operators
%  H=Operator.kron(d,n1,A1,n2,A2,...) returns the kronecker product of operators
%    A1, A2,.... The operators are positioned at subspace index n1, n2, ... in
%    the kronecker-product space. The dimension of each subspace is listed in
%    the integer vector d.
%
% Requires package:
%  - Common_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: dsum.
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 15/12/2015
% Last modified: 23/06/2017

assert(isintegervector(d) && all(d>0),...
	'QuantMech:Operator:kron:InvalidInput',...
	'Input to the subspace dimensions must be a positive integer vector.');
if nargin<=1
	H=zeros(prod(d));
	return
end
N=(nargin-1)/2;
assert(floor(N)>=N,...
	'QuantMech:Operator:kron:WrongNargin',...
	'Incorrect number of input arguments.');
D=numel(d);
j=2*(1:N);
assert(all(cellfun(@(x)isintegerscalar(x) && x>0 && x<=D,varargin(j-1))),...
	'QuantMech:Operator:kron:InvalidInput',...
	['Input to the subspace indices must be an integer scalar ',...
		'between 1 and the total number of subspaces.']);
n=[varargin{j-1}];
varargin=varargin(j);
for j=1:N
	assert(iscomplexmatrix(varargin{j}) && all(size(varargin{j})==d(n(j))),...
		'QuantMech:Operator:kron:InvalidInput',...
		['Input to the operators must be a complex square matrix whose ',...
			'dimension equals to the corresponding subspace dimension.']);
end

% Multiply operators in the same subspace
[n,j]=sort(n);
varargin=varargin(j);
j=1;
while j<N
	j=j+1;
	if n(j)<=n(j-1)
		varargin{j-1}=varargin{j-1}*varargin{j};
		n(j)=[];
		varargin(j)=[];
		j=j-1;
		N=N-1;
	end
end

% Kronecker product
H=kron(eye(prod(d(1:(n(1)-1)))),varargin{1});
j=1;
while j<N
	j=j+1;
	H=kron(H,eye(prod(d((n(j-1)+1):(n(j)-1)))));
	H=kron(H,varargin{j});
end
if n(j)<D
	H=kron(H,eye(prod(d((n(j)+1):D))));
end

end