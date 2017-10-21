function H=dsum(d,varargin)
%% Direct Sum of Operators
%  H=Operator.dsum(d,n1,A1,n2,A2,...) returns the direct sum of operators A1,
%    A2, .... The operators are positioned at subspace index n1, n2, ... in the
%    direct-sum space. The dimension of each subspace is listed in the integer
%    vector d.
%
% Requires package:
%  - Common_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: kron.
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 23/06/2017
% Last modified: 23/06/2017

assert(isintegervector(d) && all(d>0),...
	'QuantMech:Operator:dsum:InvalidInput',...
	'Input to the subspace dimensions must be a positive integer vector.');
if nargin<=1
	H=zeros(sum(d));
	return
end
N=(nargin-1)/2;
assert(floor(N)>=N,...
	'QuantMech:Operator:dsum:WrongNargin',...
	'Incorrect number of input arguments.');
D=numel(d);
i=1:N;
j=2*i;
assert(all(cellfun(@(x)isintegerscalar(x) && x>0 && x<=D,varargin(j-1))),...
	'QuantMech:Operator:dsum:InvalidInput',...
	['Input to the subspace indices must be an integer scalar ',...
		'between 1 and the total number of subspaces.']);
n=[varargin{j-1}];
varargin=varargin(j);
for j=i
	assert(iscomplexmatrix(varargin{j}) && all(size(varargin{j})==d(n(j))),...
		'QuantMech:Operator:dsum:InvalidInput',...
		['Input to the operators must be a complex square matrix whose ',...
			'dimension equals to the corresponding subspace dimension.']);
end

% Add operators in the same subspace
[n,j]=sort(n);
varargin=varargin(j);
j=1;
while j<N
	j=j+1;
	if n(j)<=n(j-1)
		varargin{j-1}=varargin{j-1}+varargin{j};
		n(j)=[];
		varargin(j)=[];
		j=j-1;
		N=N-1;
	end
end

% Direct sum
H=blkdiag(zeros(sum(d(1:(n(1)-1)))),varargin{1});
j=1;
while j<N
	j=j+1;
	H=blkdiag(H,zeros(sum(d((n(j-1)+1):(n(j)-1)))),varargin{j});
end
if n(j)<D
	H=blkdiag(H,zeros(sum(d((n(j)+1):D))));
end

end