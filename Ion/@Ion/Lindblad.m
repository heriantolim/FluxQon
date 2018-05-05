function L=Lindblad(obj,varargin)
%% Ion Lindblad Superoperator
%  L=obj.Lindlbad() returns the Lindblad superoperator that accounts for the
%    Hamiltonian and the collapse operators of the object. The basis of the
%    density matrices to which the superoperator acts on is the energy
%    eigenstates, if NumIons=1; it is the Fock states, otherwise.
%
%  L=obj.Lindblad(basis) returns the Lindblad superoperator where the density
%    matrices are in the specified basis. The options are:
%     'S' : spin basis,
%     'E' : energy eigenbasis,
%     'F' : Fock basis.
%    The 'S' and 'E' bases are applicable only when NumIons=1; and the 'F' basis
%    is applicable only when NumIons>1. Specifying 'S' or 'E', when NumIons>1,
%    will return the superperator for which NumIons=1 that acts on the density
%    matrices in the specified basis. Specifying 'F', when NumIons=1, will
%    return the superperator where the density matrices are in the 'E' basis.
%
%  L=obj.Lindblad(d,n)
%  L=obj.Lindblad(basis,d,n)
%    returns the Lindblad superoperator where the Hilbert space of the object is
%    a subspace of a larger Hilbert space. The dimensions of all the subspaces
%    is listed in the integer vector d. The subspace index for the object is
%    given by the integer scalar n.
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
% Last modified: 03/07/2017

k=1;
if k<nargin && ischar(varargin{k})
	if isscalar(varargin{k}) && any(strcmpi(varargin{k},{'S','E','F'}))
		b=upper(varargin{k});
	else
		error('FluxQon:Ion:Lindblad:InvalidInput',...
			['Input to the basis specifier must be a single character of ',...
				'either ''S'', ''E'', or ''F''']);
	end
	k=k+1;
else
	b='F';
end
if k+1<nargin
	d=varargin{k};
	assert(isintegervector(d) && numel(d)>0 && all(d>0),...
		'FluxQon:Ion:Lindblad:InvalidInput',...
		'Input to the subspace dimensions must be a positive integer vector.');
	n=varargin{k+1};
	assert(isintegerscalar(n) && n>0 && n<=numel(d),...
		'FluxQon:Ion:Lindblad:InvalidInput',...
		['Input to the subspace index must be a positive integer scalar ',...
			'that is less or equal to the number of subspaces.']);
	k=k+2;
else
	d=obj.HilbertDimension;
	n=1;
end
if k<nargin
	error('FluxQon:Ion:Lindblad:WrongNargin',...
		'Incorrect number of input arguments.');
end

D=prod(d);
L=Operator.kron(d,n,obj.Hamiltonian(b)/Constant.ReducedPlanck);
L=1i*(kron(L',eye(D))-kron(eye(D),L));

y=obj.DecayRate;
if all(y(:)==0)
	return
end

M=numel(obj.Multiplet);
MD=(2*obj.ElectronSpin+1)*(2*obj.NuclearSpin+1);
LD=MD*M;
if b=='F' && obj.NumIons>1
	if isvector(y)
		K=numel(y);
		y=reshape(y,1,K);
	else
		y=y(1,:);
		K=numel(y);
	end
	if K==1
		y=y*ones(1,LD);
	elseif K==M
		y=kron(y,ones(1,MD));
	elseif K~=LD
		error('FluxQon:Ion:Lindblad:InvalidCase',...
			['The vector dimension of the decay rate is expected to be ',...
				'equal to the number of multiplets or energy levels.']);
	end

	K=obj.Temperature;
	if K==0
		for Lj=2:LD
			C=Operator.kron(d,n,obj.Annihilation(Lj));
			L=L+y(Lj)*(kron(conj(C),C) ...
				-(kron(C.'*conj(C),eye(D))+kron(eye(D),C'*C))/2);
		end
	else
		N=obj.Energy;
		N=boltzpdf(N-N(1),K);
		N=y.*N./(1-N);
		for Lj=2:LD
			C=Operator.kron(d,n,obj.Annihilation(Lj));
			L=L+(N(Lj)+y(Lj))*(kron(conj(C),C) ...
				-(kron(C.'*conj(C),eye(D))+kron(eye(D),C'*C))/2);
			C=Operator.kron(d,n,obj.Creation(Lj));
			L=L+N(Lj)*(kron(conj(C),C) ...
				-(kron(C.'*conj(C),eye(D))+kron(eye(D),C'*C))/2);
		end
	end
else
	if isvector(y)
		y=vec2symm(y);
	end
	K=size(y,1);
	if K==1
		y=y*ones(LD);
	elseif K==M
		y=kron(y,ones(MD));
	elseif K~=LD
		error('FluxQon:Ion:Lindblad:InvalidCase',...
			['The matrix dimension of the decay rate is expected to be ',...
				'equal to the number of multiplets or energy levels.']);
	end

	K=obj.Temperature;
	if K==0
		P=1;
		LD1=1;
	else
		P=obj.Energy;
		P=boltzpdf(P-P(1),K);
		P=P/sum(P);
		LD1=LD;
	end

	if b=='S'
		W=obj.Eigenstate;
		for Li=1:LD1
			for Lj=1:LD
				if Li==Lj
					continue
				end
				C=Operator.kron(d,n,W(:,Li)*W(:,Lj)');
				L=L+P(Li)*y(Li,Lj)*(kron(conj(C),C) ...
					-(kron(C.'*conj(C),eye(D))+kron(eye(D),C'*C))/2);
			end
		end
	else
		for Li=1:LD1
			for Lj=1:LD
				if Li==Lj
					continue
				end
				C=zeros(LD);
				C(Li,Lj)=1;
				C=Operator.kron(d,n,C);
				L=L+P(Li)*y(Li,Lj)*(kron(conj(C),C) ...
					-(kron(C.'*conj(C),eye(D))+kron(eye(D),C'*C))/2);
			end
		end
	end
end

end
