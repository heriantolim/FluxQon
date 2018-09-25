function T=Transition(obj,varargin)
%% Ion Transition Operator
%  T=obj.Transition() returns a HD x HD x LD x LD matrix representing the
%    transition operators of the occupation number at the excited levels.
%    T(:,:,Li,Lj) is the operator for the transition from the level Lj to the
%    level Li, where Li,Lj=1,2,...,LD. This method is defined for calculations
%    where NumIons>1. When NumIons=1, T(;,:,Li,Lj) is a square matrix of size LD
%    x LD whose elements are all zero except that located at row Li, column Lj.
%
%  T=obj.Transition(LV1,LV2) returns a HD x HD x L1 x L2 matrix representing the
%    operators for the transitions from the levels indexed in the vector LV2 to
%    the levels indexed in the vector LV1. L1 is the length of LV1; and L2 is
%    the length of LV2.
%
%  T=obj.Transition(JV1,MV1,JV2,MV2) returns the operators for the transitions
%    from the second specified levels (JV2) to the first specified levels (JV1).
%    JV1 and JV2 are levels indexed by integers 1,2,...,MD within the respective
%    multiplets (MV1 and MV2). The MV1 and MV2 multiplets are indexed by
%    integers 1,2,...,M.
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
% See also: Annihilation, Creation, Transition.
%
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 23/06/2017
% Last modified: 23/06/2017

MD=(2*obj.ElectronSpin+1)*(2*obj.NuclearSpin+1);
LD=MD*numel(obj.Multiplet);
switch nargin
	case 1
		LV1=1:LD;
		LV2=1:LD;
	case 3
		LV1=varargin{1};
		LV2=varargin{2};
		assert(isintegervector(LV1) && all(LV1>0 & LV1<=LD) ...
				&& isintegervector(LV2) && all(LV2>0 & LV2<=LD),...
			'FluxQon:Ion:Transition:InvalidInput',...
			['Input to the level indices must be a vector whose elements are ',...
				'integers between 1 and %d.'],LD);
	case 5
		JV1=varargin{1};
		JV2=varargin{3};
		assert(isintegervector(JV1) && all(JV1>0 & JV1<=MD) ...
				&& isintegervector(JV2) && all(JV2>0 & JV2<=MD),...
			'FluxQon:Ion:Transition:InvalidInput',...
			['Input to the level indices must be a vector whose elements are ',...
				'integers between 1 and %d.'],MD);
		M=MException('FluxQon:Ion:Transition:InvalidInput',...
			['Input to the multiplet indices must be a vector whose elements ',...
				'are integers listed in the object''s Multiplet property.']);
		if isintegervector(varargin{2})
			[f,MV1]=ismember(varargin{2},obj.Multiplet);
			if ~all(f)
				throw(M);
			end
		else
			throw(M);
		end
		if isintegervector(varargin{4})
			[f,MV2]=ismember(varargin{4},obj.Multiplet);
			if ~all(f)
				throw(M);
			end
		else
			throw(M);
		end
		J=numel(JV1);
		M=numel(MV1);
		Li=0;
		LV1=zeros(1,J*M);
		for Mi=1:M
			for Ji=1:J
				Li=Li+1;
				LV1(Li)=(MV1(Mi)-1)*MD+JV1(Ji);
			end
		end
		J=numel(JV2);
		M=numel(MV2);
		Lj=0;
		LV2=zeros(1,J*M);
		for Mi=1:M
			for Ji=1:J
				Lj=Lj+1;
				LV2(Lj)=(MV2(Mi)-1)*MD+JV2(Ji);
			end
		end
	otherwise
		error('FluxQon:Ion:Transition:WrongNargin',...
			'Incorrect number of input arguments.');
end

L1=numel(LV1);
L2=numel(LV2);
if L1==0 || L2==0
	T=zeros(0,0,L1,L2);
	return
end

[LU1,~,iu1]=unique(LV1);
[LU2,~,iu2]=unique(LV2);
L1=numel(LU1);
L2=numel(LU2);
if obj.NumIons==1
	T=zeros(LD,LD,L1,L2);
	for Li=1:L1
		for Lj=1:L2
			T(LU1(Li),LU2(Lj),Li,Lj)=1;
		end
	end
else
	FD=obj.FockDimension;
	LD=LD-1;
	LU1=LU1-1;
	LU2=LU2-1;
	HD=FD^LD;
	d=FD*ones(1,LD);
	T=zeros(HD,HD,L1,L2);
	B1=Creation(FD);
	B2=Annihilation(FD);
	B3=Number(FD);
	i=1; j=1;
	if LU1(1)==0 || LU2(1)==0
		T(:,:,1,1)=obj.NumIons*eye(HD);
		for Li=1:LD
			T(:,:,1,1)=T(:,:,1,1)-Operator.kron(d,Li,B3);
		end
		if LU1(1)==0
			if L2>1
				T(:,:,1,2)=sqrt(T(:,:,1,1));
				for Lj=3:L2
					T(:,:,1,Lj)=T(:,:,1,2)*Operator.kron(d,LU2(Lj),B2);
				end
				T(:,:,1,2)=T(:,:,1,2)*Operator.kron(d,LU2(2),B2);
			end
			if LU2(1)>0
				T(:,:,1,1)=sqrt(T(:,:,1,1))*Operator.kron(d,LU2(1),B2);
			end
			i=2;
		end
		if LU2(1)==0
			if L1>1
				T(:,:,2,1)=sqrt(T(:,:,1,1));
				for Li=3:L1
					T(:,:,Li,1)=Operator.kron(d,LU1(Li),B1)*T(:,:,2,1);
				end
				T(:,:,2,1)=Operator.kron(d,LU1(2),B1)*T(:,:,2,1);
			end
			if LU1(1)>0
				T(:,:,1,1)=Operator.kron(d,LU1(1),B1)*sqrt(T(:,:,1,1));
			end
			j=2;
		end
	end
	for Li=i:L1
		for Lj=j:L2
			if LU1(Li)==LU2(Lj)
				T(:,:,Li,Lj)=Operator.kron(d,LU1(Li),B3);
			else
				T(:,:,Li,Lj)=Operator.kron(d,LU1(Li),B1,LU2(Lj),B2);
			end
		end
	end
end
T=T(:,:,iu1,iu2);

end
