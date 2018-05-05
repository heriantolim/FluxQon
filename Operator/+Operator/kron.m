function O=kron(d,varargin)
%% Kronecker Product of Operators
%  O=Operator.kron(d,n1,A1,n2,A2,...) returns the kronecker product of the
%  operators A1, A2,.... The operators are positioned at subspace index n1, n2,
%  ... in the kronecker-product space. The dimension of each subspace must be
%  listed in the integer vector d.
%
% Requires package:
%  - MatCommon_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: dsum.
%
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 15/12/2015
% Last modified: 23/06/2017

assert(isintegervector(d) && all(d>0),...
	'FluxQon:Operator:kron:InvalidInput',...
	'Input to the subspace dimensions must be a positive integer vector.');
if nargin<=1
	O=zeros(prod(d));
	return
end
N=(nargin-1)/2;
assert(floor(N)>=N,...
	'FluxQon:Operator:kron:WrongNargin',...
	'Incorrect number of input arguments.');
D=numel(d);
j=2*(1:N);
assert(all(cellfun(@(x)isintegerscalar(x) && x>0 && x<=D,varargin(j-1))),...
	'FluxQon:Operator:kron:InvalidInput',...
	['Input to the subspace indices must be an integer scalar ',...
		'between 1 and the total number of subspaces.']);
n=[varargin{j-1}];
varargin=varargin(j);
for j=1:N
	assert(iscomplexmatrix(varargin{j}) && all(size(varargin{j})==d(n(j))),...
		'FluxQon:Operator:kron:InvalidInput',...
		['Input to the operators must be a complex square matrix whose ',...
			'dimension equals to the corresponding subspace dimension.']);
end

% Multiply the operators in the same subspace.
[n,j]=sort(n);
varargin=varargin(j);
j=1;
while j<N
	j=j+1;
	if n(j)<=n(j-1)
		varargin{j-1}=varargin{j-1}*varargin{j};% matrix multiplication
		n(j)=[];
		varargin(j)=[];
		j=j-1;
		N=N-1;
	end
end

% Kronecker product.
O=kron(eye(prod(d(1:(n(1)-1)))),varargin{1});
j=1;
while j<N
	j=j+1;
	O=kron(O,eye(prod(d((n(j-1)+1):(n(j)-1)))));
	O=kron(O,varargin{j});
end
if n(j)<D
	O=kron(O,eye(prod(d((n(j)+1):D))));
end

end