function L=Lindblad(obj,varargin)
%% Qubit Lindblad Superoperator
%  L=obj.Lindlbad() returns the Lindblad superoperator that accounts for the
%    Hamiltonian and the collapse operators of the object.
%
%  L=obj.Lindblad(d,n) returns the Lindblad superoperator where the Hilbert
%    space of the object is a subspace of a larger Hilbert space. The dimensions
%    of all the subspaces is listed in the integer vector d. The subspace index
%    for the object is given by the integer scalar n.
%
% Requires package:
%  - MatCommon_v1.0.0+
%  - PhysConst_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
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

y=obj.DecayRate;
if y>0
	T=obj.Temperature;
	if T==0
		C=Operator.kron(d,n,Pauli.M);
		L=L+y*(kron(conj(C),C)-(kron(C.'*conj(C),eye(D))+kron(eye(D),C'*C))/2);
	else
		P=[1;boltzpdf(obj.Energy,T)];
		P=P/sum(P);
		C=Operator.kron(d,n,Pauli.M);
		L=L+P(1)*y*(kron(conj(C),C) ...
			-(kron(C.'*conj(C),eye(D))+kron(eye(D),C'*C))/2);
		C=Operator.kron(d,n,Pauli.P);
		L=L+P(2)*y*(kron(conj(C),C) ...
			-(kron(C.'*conj(C),eye(D))+kron(eye(D),C'*C))/2);
	end
end

y=obj.DephasingRate;
if y>0
	C=Operator.kron(d,n,Pauli.Z);
	L=L+y*(kron(conj(C),C)-(kron(C.'*conj(C),eye(D))+kron(eye(D),C'*C))/2);
end

end
