function C=Creation(obj,varargin)
%% Ion Creation Operator
%  C=obj.Creation() returns a HD x HD x LD matrix representing the creation
%    operators of the occupation number at the excited levels. C(:,:,Li) is the
%    creation operator for the level Li, where Li=1,2,...,LD. C(:,:,1) is
%    not a creation operator for the ground level, but rather a self-adjoint
%    operator defined in Ref. [1]. This method is defined for calculations where
%    NumIons>1. When NumIons=1, C(;,:,Li) is a square matrix of size LD x LD
%    whose elements are all zero except that located at row Li, column 1.
%
%  C=obj.Creation(LV) returns the creation operators for the levels specified in
%    the input vector LV. LV is a vector whose elements are integers between 1
%    to LD. The output C has a dimension of HD x HD x L, where L is the length
%    of LV.
%
%  C=obj.Creation(JV,MV) returns the creation operators for the specified levels
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
%  - MatCommon_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: Annihilation, Number, Transition.
%
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 20/06/2017
% Last modified: 23/06/2017

MD=(2*obj.ElectronSpin+1)*(2*obj.NuclearSpin+1);
LD=MD*numel(obj.Multiplet);
switch nargin
	case 1
		LV=1:LD;
	case 2
		LV=varargin{1};
		assert(isintegervector(LV) && all(LV>0 & LV<=LD),...
			'FluxQon:Ion:Creation:InvalidInput',...
			['Input to the level indices must be a vector whose elements are ',...
				'integers between 1 and %d.'],LD);
	case 3
		JV=varargin{1};
		assert(isintegervector(JV) && all(JV>0 & JV<=MD),...
			'FluxQon:Ion:Creation:InvalidInput',...
			['Input to the level indices must be a vector whose elements are ',...
				'integers between 1 and %d.'],MD);
		M=MException('FluxQon:Ion:Creation:InvalidInput',...
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
		error('FluxQon:Ion:Creation:WrongNargin',...
			'Incorrect number of input arguments.');
end

if numel(LV)==0
	C=[];
	return
end

[LU,~,iu]=unique(LV);
L=numel(LU);
if obj.NumIons==1
	C=zeros(LD,LD,L);
	for Li=1:L
		C(LU(Li),1,Li)=1;
	end
else
	FD=obj.FockDimension;
	LD=LD-1;
	LU=LU-1;
	HD=FD^LD;
	d=FD*ones(1,LD);
	C=zeros(HD,HD,L);
	if LU(1)==0
		B=Number(FD);
		C(:,:,1)=obj.NumIons*eye(HD);
		for Li=1:LD
			C(:,:,1)=C(:,:,1)-Operator.kron(d,Li,B);
		end
		C(:,:,1)=sqrt(C(:,:,1));
		i=2;
	else
		i=1;
	end
	B=Creation(FD);
	for Li=i:L
		C(:,:,Li)=Operator.kron(d,LU(Li),B);
	end
end
C=C(:,:,iu);

end
