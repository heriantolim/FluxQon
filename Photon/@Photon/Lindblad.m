function L=Lindblad(obj,varargin)
%% Photon Lindblad Superoperator
%  L=obj.Lindlbad() returns the Lindblad superoperator that accounts for the
%    Hamiltonian and the collapse operators of the object.
%
%  L=obj.Lindblad(d,n) returns the Lindblad superoperator where the Hilbert
%    space of the object is a subspace of a larger Hilbert space. The dimensions
%    of all the subspaces is listed in the integer vector d. The subspace index
%    for the object is given by the integer scalar n.
%
% Requires package:
%  - Common_v1.0.0+
%  - PhysConst_v1.0.0+
%  - QuantMech_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 18/06/2017
% Last modified: 19/06/2017

if nargin==1
	d=obj.HilbertDimension;
	n=1;
elseif nargin==3
	d=varargin{1};
	n=varargin{2};
	assert(isintegervector(d) && numel(d)>0 && all(d>0),...
		'FluxQon:Qubit:Lindblad:InvalidInput',...
		'Input to the subspace dimensions must be a positive integer vector.');
	assert(isintegerscalar(n) && n>0 && n<=numel(d),...
		'FluxQon:Qubit:Lindblad:InvalidInput',...
		['Input to the subspace index must be a positive integer scalar ',...
			'that is less or equal to the number of subspaces.']);
else
	error('FluxQon:Qubit:Lindblad:WrongNargin',...
		'Incorrect number of input arguments.');
end

D=prod(d);
L=Operator.kron(d,n,obj.Hamiltonian/Constant.ReducedPlanck);
L=1i*(kron(L',eye(D))-kron(eye(D),L));

N=boltzpdf(obj.Energy,obj.Temperature);
N=N/(1-N);

K=obj.DecayRate;
if K>0
	C=Operator.kron(d,n,obj.Annihilation);
	L=L+K*(N+1)*(kron(conj(C),C)-(kron(C.'*conj(C),eye(D))+kron(eye(D),C'*C))/2);
end

P=obj.PumpRate;
if P>0 || K>0
	C=Operator.kron(d,n,obj.Creation);
	L=L+(P+K*N)*(kron(conj(C),C)-(kron(C.'*conj(C),eye(D))+kron(eye(D),C'*C))/2);
end

end