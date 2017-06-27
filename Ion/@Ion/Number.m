function N=Number(obj,varargin)
%% Ion Number Operator
%  N=obj.Number() returns a HD x HD x LD matrix representing the number
%    operators of the occupation number at the excited levels. N(:,:,Li) is the
%    number operator for the level Li, where Li=1,2,...,LD. This method is
%    defined for calculations where NumIons>1. When NumIons=1, N(;,:,Li) is a
%    square matrix of size LD x LD whose elements are all zero except that
%    located at row Li, column Li.
%
%  N=obj.Number(LV) returns the number operators for the levels specified in the
%    input vector LV. LV is a vector whose elements are integers between 1 to
%    LD. The output C has a dimension of HD x HD x L, where L is the length of
%    LV.
%
%  N=obj.Number(JV,MV) returns the number operators for the specified levels
%    (JV) in the specified multiplets (MV). The multiplets are indexed by
%    integers 1,2,...,M, and the levels in each multiplets are indexed by
%    integers 1,2,...,MD. JV and MV are a vector of the indices.
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
% References:
%  [1] Z. Kurucz and K. Molmer, Phys. Rev. A 81, 1 (2010).
%
% Requires package:
%  - Common_v1.0.0+
%  - QuantMech_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: Annihilation, Creation, Transition.
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 22/06/2017
% Last modified: 23/06/2017

MD=(2*obj.ElectronSpin+1)*(2*obj.NuclearSpin+1);
LD=MD*numel(obj.Multiplet);
switch nargin
	case 1
		LV=1:LD;
	case 2
		LV=varargin{1};
		assert(isintegervector(LV) && all(LV>0 & LV<=LD),...
			'FluxQon:Ion:Number:InvalidInput',...
			['Input to the level indices must be a vector whose elements are ',...
				'integers between 1 and %d.'],LD);
	case 3
		JV=varargin{1};
		assert(isintegervector(JV) && all(JV>0 & JV<=MD),...
			'FluxQon:Ion:Number:InvalidInput',...
			['Input to the level indices must be a vector whose elements are ',...
				'integers between 1 and %d.'],MD);
		M=MException('FluxQon:Ion:Number:InvalidInput',...
			['Input to the multiplet indices must be a vector whose elements ',...
				'are integers listed in the object''s Multiplet property.']);
		if isintegervector(varargin{2})
			[f,MV]=ismember(varargin{2},obj.Multiplet);
			if ~all(f)
				throw(M);
			end
		else
			throw(M);
		end
		J=numel(JV);
		M=numel(MV);
		Li=0;
		LV=zeros(1,J*M);
		for Mi=1:M
			for Ji=1:J
				Li=Li+1;
				LV(Li)=(MV(Mi)-1)*MD+JV(Ji);
			end
		end
	otherwise
		error('FluxQon:Ion:Number:WrongNargin',...
			'Incorrect number of input arguments.');
end

if numel(LV)==0
	N=[];
	return
end

[LU,~,iu]=unique(LV);
L=numel(LU);
if obj.NumIons==1
	N=zeros(LD,LD,L);
	for Li=1:L
		N(LU(Li),LU(Li),Li)=1;
	end
else
	FD=obj.FockDimension;
	LD=LD-1;
	LU=LU-1;
	HD=FD^LD;
	d=FD*ones(1,LD);
	N=zeros(HD,HD,L);
	B=Number(FD);
	if LU(1)==0
		N(:,:,1)=obj.NumIons*eye(HD);
		for Li=1:LD
			N(:,:,1)=N(:,:,1)-Operator.kron(d,Li,B);
		end
		i=2;
	else
		i=1;
	end
	for Li=i:L
		N(:,:,Li)=Operator.kron(d,LU(Li),B);
	end
end
N=N(:,:,iu);

end