function H=HyperfineHamiltonian(obj,varargin)
%% Ion Hyperfine Hamiltonian
%  H=obj.HyperfineHamiltonian() returns the hyperfine Hamiltonian in the basis
%    of the energy eigenstates, if NumIons=1, or in the basis of the Fock
%    states, otherwise.
%
%  H=obj.ZeemanHamiltonian(basis) returns the hyperfine Hamiltonian in the
%    specified basis. The options are:
%     'S' : spin basis,
%     'E' : energy eigenbasis,
%     'F' : Fock basis.
%    The 'S' and 'E' bases are applicable only when NumIons=1; and the 'F' basis
%    is applicable only when NumIons>1. Specifying 'S' or 'E', when NumIons>1,
%    will return the Hamiltonian for which NumIons=1 in the specified basis.
%    Specifying 'F', when NumIons=1, will return the Hamiltonian in the 'E'
%    basis.
%
% Requires package:
%  - Common_v1.0.0+
%  - QuantMech_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: ZeemanHamiltonian, QuadrupoleHamiltonian, Hamiltonian.
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 17/06/2017
% Last modified: 23/06/2017

switch nargin
	case 1
		b='F';
	case 2
		if ischar(varargin{1}) && isscalar(varargin{1}) ...
				&& any(strcmpi(varargin{1},{'S','E','F'}))
			b=upper(varargin{1});
		else
			error('FluxQon:Ion:HyperfineHamiltonian:InvalidInput',...
				['Input to the basis specifier must be a single character of ',...
					'either ''S'', ''E'', or ''F''']);
		end
	otherwise
		error('FluxQon:Ion:HyperfineHamiltonian:WrongNargin',...
			'Incorrect number of input arguments.');
end

MV=obj.Multiplet;
M=numel(MV);
LD=obj.NumLevels;
T=obj.HyperfineTensor;

if all(T(:)==0)
	if b=='F' && obj.NumIons>1
		H=zeros(obj.FockDimension^(LD-1));
	else
		H=zeros(LD);
	end
	return
end

H=zeros(LD);
f=size(T,3)==1;
for k=1:M
	if f
		Ti=T;
	else
		Ti=T(:,:,k);
	end
	if any(Ti(:)~=0)
		if isscalar(Ti)
			s={'X','Y','Z'};
			for i=1:3
				H=H+Ti*obj.ElectronSpinOperator('S',s{i},MV(k)) ...
					*obj.NuclearSpinOperator('S',s{i},MV(k));
			end
		else
			A=obj.ElectronSpinOperator('S',MV(k));
			C=obj.NuclearSpinOperator('S',MV(k));
			for i=1:3
				for j=1:3
					H=H+Ti(i,j)*A(:,:,i)*C(:,:,j);
				end
			end
		end
	end
end

if b~='S'
	A=obj.Eigenstate;
	H=A'*H*A;
end

if b=='F' && obj.NumIons>1
	B=H-obj.Energy(1)*eye(LD);
	if B(1,1)==0
		H=zeros(obj.FockDimension^(LD-1));
	else
		H=-obj.NumIons*B(1,1)*eye(obj.FockDimension^(LD-1));
	end
	A=obj.Annihilation;
	C=obj.Creation;
	for i=1:LD
		for j=1:LD
			H=H+C(:,:,i)*B(i,j)*A(:,:,j);
		end
	end
end

b=obj.RoundRelTol;
if b>0
	try
		H=round2herm(H,b);
	catch ME1
		if isempty(regexpi(ME1.identifier,':ToleranceExceeded$','once'))
			rethrow(ME1);
		else
			warning(...
				'FluxQon:Ion:HyperfineHamiltonian:NonHermitian',...
				'The hyperfine Hamiltonian is non Hermitian.');
		end
	end
end

end
