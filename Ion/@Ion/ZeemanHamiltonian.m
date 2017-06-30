function H=ZeemanHamiltonian(obj,varargin)
%% Ion Zeeman Hamiltonian
%  H=obj.ZeemanHamiltonian() returns the Zeeman Hamiltonian in the basis of the
%    energy eigenstates, if NumIons=1, or in the basis of the Fock states,
%    otherwise.
%
%  H=obj.ZeemanHamiltonian(B) returns the Zeeman Hamiltonian where B is used as
%    the magnetic field instead of the object's MagneticField property. B can be
%    specified as either a real scalar or a real vector of length 3. When
%    specified as a real scalar, B is treated as being equal to [0;0;B].
%
%  H=obj.ZeemanHamiltonian(basis)
%  H=obj.ZeemanHamiltonian(basis,B) returns the Zeeman Hamiltonian in the
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
%  - PhysConst_v1.0.0+
%  - QuantMech_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: HyperfineHamiltonian, QuadrupoleHamiltonian, Hamiltonian.
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 17/06/2017
% Last modified: 23/06/2017

k=1;
if k<nargin && ischar(varargin{k})
	if isscalar(varargin{k}) && any(strcmpi(varargin{k},{'S','E','F'}))
		b=upper(varargin{k});
	else
		error('FluxQon:Ion:ZeemanHamiltonian:InvalidInput',...
			['Input to the basis specifier must be a single character of ',...
				'either ''S'', ''E'', or ''F''']);
	end
	k=k+1;
else
	b='F';
end
if k<nargin
	if isrealscalar(varargin{k})
		B=[0;0;varargin{k}];
	elseif isrealvector(varargin{k}) && numel(varargin{k})==3
		B=reshape(varargin{k},3,1);
	else
		error('FluxQon:Ion:ZeemanHamiltonian:InvalidInput',...
			['Input to the magnetic field must be either ',...
				'a real scalar or a real vector of length 3.']);
	end
	k=k+1;
else
	B=obj.MagneticField;
end
if k<nargin
	error('FluxQon:Ion:ZeemanHamiltonian:WrongNargin',...
		'Incorrect number of input arguments.');
end

MV=obj.Multiplet;
M=numel(MV);
LD=obj.NumLevels;
T1=Constant.BohrMagneton*obj.ElectronZeemanTensor;
T2=Constant.NuclearMagneton*obj.NuclearZeemanTensor;

f1=all(T1(:)==0);
f2=all(T2(:)==0);
if f1 && f2
	if b=='F' && obj.NumIons>1
		H=zeros(obj.FockDimension^(LD-1));
	else
		H=zeros(LD);
	end
	return
end

% Electron Zeeman
H=zeros(LD);
if ~f1
	f1=size(T1,3)==1;
	for k=1:M
		if f1
			Ti=T1;
		else
			Ti=T1(:,:,k);
		end
		if any(Ti(:)~=0)
			if isscalar(Ti)
				s={'X','Y','Z'};
				for i=1:3
					if B(i)~=0
						H=H+Ti*B(i)*obj.ElectronSpinOperator('S',s{i},MV(k));
					end
				end
			else
				A=obj.ElectronSpinOperator('S',MV(k));
				for i=1:3
					if B(i)~=0
						for j=1:3
							H=H+B(i)*Ti(i,j)*A(:,:,j);
						end
					end
				end
			end
		end
	end
end

% Nuclear Zeeman
if ~f2
	f2=size(T2,3)==1;
	for k=1:M
		if f2
			Ti=T2;
		else
			Ti=T2(:,:,k);
		end
		if any(Ti(:)~=0)
			if isscalar(Ti)
				s={'X','Y','Z'};
				for i=1:3
					if B(i)~=0
						H=H-Ti*B(i)*obj.NuclearSpinOperator('S',s{i},MV(k));
					end
				end
			else
				A=obj.NuclearSpinOperator('S',MV(k));
				for i=1:3
					if B(i)~=0
						for j=1:3
							H=H-B(i)*Ti(i,j)*A(:,:,j);
						end
					end
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
				'FluxQon:Ion:ElectronZeemanHamiltonian:NonHermitian',...
				'The Zeeman Hamiltonian is non Hermitian.');
		end
	end
end

end