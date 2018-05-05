function L=Lindblad(varargin)
%% Lindblad Superoperator
%  This function returns the superoperator L for the Lindblad master equation,
%  which is generally characterized by a Hamiltonian H and a set of collapse
%  operators Ci with a collapse rate ri.
%
%  L=Lindblad(H) returns the superoperator defined by only the Hamiltonian H.
%
%  L=Lindblad(r1,C1,r2,C2,...) returns the superoperator defined by the collapse
%    operators Ci with collapse rate ri. All the Ci must be a square complex
%    matrix and have the same dimensions.
%
%  L=Lindblad(H,r1,C1,r2,C2,...) returns the superoperator that is defined by
%    both the Hamiltonian H and the collapse operators Ci with collapse rate ri.
%    Each Ci must have the same dimensions as H.
%
%  L=Lindblad(d,n1,r1,C1,n2,r2,C2,...) returns the superoperator where the
%    kronecker product dimension is specified by the vector d, and each collapse
%    operator Ci is positioned at index ni in d and associated with a collapse
%    rate ri.
%
%  L=Lindblad(H,d,n1,r1,C1,n2,r2,C2,...) returns the superoperator defined by
%    the Hamiltonian H and the collapse operators Ci.
%
% Inputs:
%  - H  : Hamiltonian, specified as a Hermitian complex matrix.
%  - d  : Subspace dimensions, specified as an integer vector.
%  - ni : Subspace indices, specified as a positive integer scalar.
%  - ri : Collapse rate, specified as a positive real scalar.
%  - Ci : Collapse operator, specified as a square complex matrix.
%
% Outputs:
%  - L  : Lindblad superoperator, returned as a square complex matrix. If the
%         dimension of the Hamiltonian is D, then L has a dimension of D^2.
%
% Requires package:
%  - MatCommon_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: TDSE.
%
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 15/06/2017
% Last modified: 15/06/2017

K=nargin;
assert(K>0,...
	'FluxQon:Operator:Lindblad:WrongNargin',...
	'Incorrect number of input arguments.');

k=1;
if iscomplexmatrix(varargin{k}) && all(size(varargin{k})>1)
	assert(ishermitian(varargin{k}),...
		'FluxQon:Operator:Lindblad:InvalidInput',...
		'Input to the Hamiltonian matrix must be Hermitian.');
	L=varargin{k}/Constant.ReducedPlanck;
	D=size(L,1);
	L=1i*(kron(L',eye(D))-kron(eye(D),L));
	k=k+1;
else
	L=0;
	D=0;
end

if k+1>K
	return
end

if numel(varargin{k})==1
	flag=false;
else
	flag=true;
	assert(isintegervector(varargin{k}) && all(varargin{k}>0),...
		'FluxQon:Operator:Lindblad:InvalidInput',...
		'Input to the subspace dimensions must be a positive integer vector.');
	d=varargin{k};
	if D==0
		D=prod(d);
	elseif D~=prod(d)
		error('FluxQon:Operator:Lindblad:InvalidInput',...
			['The dimension of the Hamiltonian must equal to that of the ',...
				'kronecker product space of the collapse operators.']);
	end
	k=k+1;
end

while k<K
	if flag
		n=varargin{k};
		k=k+1;
	end

	assert(isrealscalar(varargin{k}) && varargin{k}>0,...
		'FluxQon:Operator:Lindblad:InvalidInput',...
		'Input to the collapse rate must be a positive real scalar.');
	r=varargin{k};
	k=k+1;

	if flag
		try
			C=Operator.kron(d,n,varargin{k});
		catch ME1
			if isempty(regexpi(ME1.identifier,':InvalidInput$','once'))
				rethrow(ME1);
			else
				ME=MException('FluxQon:Operator:Lindblad:InvalidInput',...
					'Input argument %d and %d failed validation.',k-2,k);
				ME=addCause(ME,ME1);
				throw(ME);
			end
		end
	else
		C=size(varargin{k});
		if D==0
			D=C(1);
		end
		assert(iscomplexmatrix(varargin{k}) && all(C>0) && all(C==D),...
			'FluxQon:Operator:Lindblad:InvalidInput',...
			['Each collapse operator must be a square complex matrix and ',...
				'have the same dimension (as the Hamiltonian).']);
		C=varargin{k};
	end
	k=k+1;

	L=L+r*(kron(conj(C),C)-(kron(C.'*conj(C),eye(D))+kron(eye(D),C'*C))/2);
end

if k<=K
	warning('FluxQon:Operator:Lindblad:IgnoredInput',...
		'An extra input argument was provided, but is ignored.');
end

end
