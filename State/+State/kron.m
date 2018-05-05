function psi=kron(d,varargin)
%% Kronecker Product of States
%  psi=State.kron(d,n1,psi1,n2,psi2,...) returns the kronecker product of the
%  states psi1, psi2, ... The states will be normalized before being combined in
%  the kronecker product. The states are positioned at subspace index n1, n2,
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
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 21/05/2017
% Last modified: 23/06/2017

assert(isintegervector(d) && all(d>0),...
	'FluxQon:State:kron:InvalidInput',...
	'Input to the subspace dimensions must be a positive integer vector.');
if nargin<=1
	psi=zeros(prod(d),1);
	return
end
N=(nargin-1)/2;
assert(floor(N)>=N,...
	'FluxQon:State:kron:WrongNargin',...
	'Incorrect number of input arguments.');
D=numel(d);
j=2*(1:N);
assert(all(cellfun(@(x)isintegerscalar(x) && x>0 && x<=D,varargin(j-1))),...
	'FluxQon:State:kron:InvalidInput',...
	['Input to the subspace indices must be an integer scalar ',...
		'between 1 and the total number of subspaces.']);
n=[varargin{j-1}];
varargin=varargin(j);
for j=1:N
	assert(iscomplexvector(varargin{j}) && numel(varargin{j})==d(n(j)),...
		'FluxQon:State:kron:InvalidInput',...
		['Input to the states must be a complex vector whose ',...
			'length equals to the corresponding subspace dimension.']);
		
	% Normalize the states.
	varargin{j}=normc(varargin{j}(:));
end

% Multiply the states in the same subspace.
[n,j]=sort(n);
varargin=varargin(j);
j=1;
while j<N
	j=j+1;
	if n(j)<=n(j-1)
		varargin{j-1}=varargin{j-1}.*varargin{j};
		n(j)=[];
		varargin(j)=[];
		j=j-1;
		N=N-1;
	end
end

% Kronecker Product.
psi=kron(eye(prod(d(1:(n(1)-1))),1),varargin{1});
j=1;
while j<N
	j=j+1;
	psi=kron(psi,eye(prod(d((n(j-1)+1):(n(j)-1))),1));
	psi=kron(psi,varargin{j});
end
if n(j)<D
	psi=kron(psi,eye(prod(d((n(j)+1):D)),1));
end

end