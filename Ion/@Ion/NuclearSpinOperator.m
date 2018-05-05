function IO=NuclearSpinOperator(obj,varargin)
%% Ion Nuclear Spin Operator
%  IO=obj.NuclearSpinOperator() returns a HD x HD x 3 x M matrix. Each
%    submatrice IO(:,:,1,Mi), IO(:,:,2,Mi), and IO(:,:,3,Mi) represent the
%    nuclear spin-X, -Y, and -Z operators, respectively, for the multiplet
%    obj.Multiplet(Mi).
%
%  IO=obj.NuclearSpinOperator(basis) returns the nuclear spin-X, -Y, and -Z
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
%  IO=obj.NuclearSpinOperator('Y','X',...) returns nuclear the spin-Y, -X, ...
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
%  IO=obj.NuclearSpinOperator(basis,'Y','X',...) returns the nuclear spin-Y, -X,
%    ... in the specified basis.
%
%  IO=obj.NuclearSpinOperator(...,MV) returns the nuclear spin operators only
%    for the multiplets specified in the integer vector MV.
%
% Glossary:
%  - S : Nuclear spin.
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
% See also: ElectronSpinOperator.
%
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 17/06/2017
% Last modified: 24/06/2017

b='S';
K=nargin-1;
if K>0 && ~ischar(varargin{K})
	M=MException('FluxQon:Ion:NuclearSpinOperator:InvalidInput',...
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
		'FluxQon:Ion:NuclearSpinOperator:InvalidInput',...
		'All but the last input arguments must be single characters.');
	if any(strcmpi(varargin{1},{'S','E','F'}))
		b=upper(varargin{1});
		varargin(1)=[];
		K=K-1;
	end
end
if K>0
	assert(all(cellfun(@(x)any(strcmp(x,{'X','Y','Z','P','M'})),varargin)),...
		'FluxQon:Ion:NuclearSpinOperator:InvalidInput',...
		['Specifiers for the type of the spin operators must be either ',...
			'''X'', ''Y'', ''Z'', ''P'', or ''M''.']);
else
	K=3;
	varargin={'X','Y','Z'};
end

I=obj.NuclearSpin;
SD=2*obj.ElectronSpin+1;
MD=SD*(2*I+1);
d=numel(obj.Multiplet);
LD=MD*d;

if I==0
	if b=='F' && obj.NumberIons>1
		HD=obj.FockDimension^(LD-1);
		IO=zeros(HD,HD,K,M);
	else
		IO=zeros(LD,LD,K,M);
	end
	return
end

d=MD*ones(1,d);
IO=zeros(LD,LD,K,M);
for Mi=1:M
	for k=1:K
		IO(:,:,k,Mi)=Operator.dsum(d,MV(Mi),kron(eye(SD),Spin.(varargin{k})(I)));
	end
end

if b~='S'
	W=obj.Eigenstate;
	for Mi=1:M
		for k=1:K
			IO(:,:,k,Mi)=W'*IO(:,:,k,Mi)*W;
		end
	end
end

if b=='F' && obj.NumIons>1
	HD=obj.FockDimension^(LD-1);
	I=IO;
	if I(1,1)==0
		IO=zeros(HD,HD,K,M);
	else
		IO=repmat(-obj.NumIons*I(1,1)*eye(HD),1,1,K,M);
	end
	A=obj.Annihilation;
	C=obj.Creation;
	for Mi=1:M
		for k=1:K
			for Li=1:LD
				for Lj=1:LD
					IO(:,:,k,Mi)=IO(:,:,k,Mi)+C(:,:,Li)*I(Li,Lj)*A(:,:,Lj);
				end
			end
		end
	end
end

end
