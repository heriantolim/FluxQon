function SO=ElectronSpinOperator(obj,varargin)
%% Ion Electron Spin Operator
%  SO=obj.ElectronSpinOperator() returns a HD x HD x 3 x M matrix. Each
%    submatrice SO(:,:,1,Mi), SO(:,:,2,Mi), and SO(:,:,3,Mi) represent the
%    electron spin-X, -Y, and -Z operators, respectively, for the multiplet
%    obj.Multiplet(Mi).
%
%  SO=obj.ElectronSpinOperator(basis) returns the electron spin-X, -Y, and -Z
%    operators in the specified basis. The options are:
%     'S' : spin basis,
%     'E' : energy eigenbasis,
%     'F' : Fock basis.
%    The 'S' and 'E' bases are applicable only when NumIons=1; and the 'F' basis
%    is applicable only when NumIons>1. Specifying 'S' or 'E', when NumIons>1,
%    will return the spin operators for which NumIons=1 in the specified basis.
%    Specifying 'F', when NumIons=1, will return the spin operators in the 'E'
%    basis.
%
%  SO=obj.ElectronSpinOperator('Y','X',...) returns the electron spin-Y, -X, ...
%    operators concatenated along the third dimension. The basis of the
%    operators is either the energy eigenstates or the Fock states, depending on
%    the value of NumIons. The input arguments are single characters specifying
%    the type of the spin operators to be returned. The options are:
%     'X' : Spin-X,
%     'Y' : Spin-Y,
%     'Z' : Spin-Z,
%     'P' : Spin raising,
%     'M' : Spin lowering.
%    The order of the concatenation matches the order of the specified types.
%
%  SO=obj.ElectronSpinOperator(basis,'Y','X',...) returns the electron spin-Y,
%    -X, ... in the specified basis.
%
%  SO=obj.ElectronSpinOperator(...,MV) returns the electron spin operators only
%    for the multiplets specified in the integer vector MV.
%
% Glossary:
%  - S : Electron spin.
%  - I : Nuclear spin.
%  - M : Number of multiplets.
%  - SD=2*S+1 : Dimension of each electron-spin subspace.
%  - ID=2*I+1 : Dimension of each nuclear-spin subspace.
%  - MD=SD*ID : Dimension of each multiplet.
%  - LD=MD*M  : Total dimension of all multiplets.
%  - FD       : Dimension of each Fock subspace.
%  - HD<=FD^(LD-1) : Dimension of the Hilbert space.
%
% Requires package:
%  - MatCommon_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: NuclearSpinOperator.
%
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 17/06/2017
% Last modified: 24/06/2017

b='S';
K=nargin-1;
if K>0 && ~ischar(varargin{K})
	M=MException('FluxQon:Ion:ElectronSpinOperator:InvalidInput',...
		['Input to the multiplet indices must be a vector whose elements ',...
			'are integers listed in the object''s Multiplet property.']);
	if isintegervector(varargin{K})
		[f,MV]=ismember(varargin{K},obj.Multiplet);
		if ~all(f)
			throw(M);
		end
	else
		throw(M);
	end
	M=numel(MV);
	varargin(K)=[];
	K=K-1;
else
	M=numel(obj.Multiplet);
	MV=1:M;
end
if K>0
	assert(all(cellfun(@(x)ischar(x) && isscalar(x),varargin)),...
		'FluxQon:Ion:ElectronSpinOperator:InvalidInput',...
		'All but the last input arguments must be single characters.');
	if any(strcmpi(varargin{1},{'S','E','F'}))
		b=upper(varargin{1});
		varargin(1)=[];
		K=K-1;
	end
end
if K>0
	assert(all(cellfun(@(x)any(strcmp(x,{'X','Y','Z','P','M'})),varargin)),...
		'FluxQon:Ion:ElectronSpinOperator:InvalidInput',...
		['Specifiers for the type of the spin operators must be either ',...
			'''X'', ''Y'', ''Z'', ''P'', or ''M''.']);
else
	K=3;
	varargin={'X','Y','Z'};
end

S=obj.ElectronSpin;
ID=2*obj.NuclearSpin+1;
MD=(2*S+1)*ID;
d=numel(obj.Multiplet);
LD=MD*d;

if S==0
	if b=='F' && obj.NumberIons>1
		HD=obj.FockDimension^(LD-1);
		SO=zeros(HD,HD,K,M);
	else
		SO=zeros(LD,LD,K,M);
	end
	return
end

d=MD*ones(1,d);
SO=zeros(LD,LD,K,M);
for Mi=1:M
	for k=1:K
		SO(:,:,k,Mi)=Operator.dsum(d,MV(Mi),kron(Spin.(varargin{k})(S),eye(ID)));
	end
end

if b~='S'
	W=obj.Eigenstate;
	for Mi=1:M
		for k=1:K
			SO(:,:,k,Mi)=W'*SO(:,:,k,Mi)*W;
		end
	end
end

if b=='F' && obj.NumIons>1
	HD=obj.FockDimension^(LD-1);
	S=SO;
	if S(1,1)==0
		SO=zeros(HD,HD,K,M);
	else
		SO=repmat(-obj.NumIons*S(1,1)*eye(HD),1,1,K,M);
	end
	A=obj.Annihilation;
	C=obj.Creation;
	for Mi=1:M
		for k=1:K
			for Li=1:LD
				for Lj=1:LD
					SO(:,:,k,Mi)=SO(:,:,k,Mi)+C(:,:,Li)*S(Li,Lj)*A(:,:,Lj);
				end
			end
		end
	end
end

end
