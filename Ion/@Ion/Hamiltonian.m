function H=Hamiltonian(obj,varargin)
%% Ion Hamiltonian
%  Calling this method resets the object's Energy and Eigenstate property.
%
%  H=obj.Hamiltonian() returns the Hamiltonian in the basis of the energy
%    eigenstates, if NumIons=1, or in the basis of the Fock states, otherwise.
%
%  H=obj.Hamiltonian(basis) returns the Hamiltonian in the specified basis.
%    The options are:
%     'S' : spin basis,
%     'E' : energy eigenbasis,
%     'F' : Fock basis.
%    The 'S' and 'E' bases are applicable only when NumIons=1; and the 'F' basis
%    is applicable only when NumIons>1. Specifying 'S' or 'E', when NumIons>1,
%    will return the Hamiltonian for which NumIons=1 in the specified basis.
%    Specifying 'F', when NumIons=1, will return the Hamiltonian in the 'E'
%    basis.
%
%  H=obj.Hamiltonian(B)
%  H=obj.Hamiltonian(basis,B)
%    returns the Hamiltonian where B is used as the magnetic field instead of
%    the object's MagneticField property. B can be specified as either a real
%    scalar or a real vector of length 3. When specified as a real scalar, B is
%    treated as being equal to [0;0;B].
%
% Requires package:
%  - MatCommon_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: Hamiltonian, HyperfineHamiltonian, QuadrupoleHamiltonian.
%
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 17/06/2017
% Last modified: 24/06/2017

k=1;
if k<nargin && ischar(varargin{k})
	if isscalar(varargin{k}) && any(strcmpi(varargin{k},{'S','E','F'}))
		b=upper(varargin{k});
	else
		error('FluxQon:Ion:Hamiltonian:InvalidInput',...
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
		error('FluxQon:Ion:Hamiltonian:InvalidInput',...
			['Input to the magnetic field must be either ',...
				'a real scalar or a real vector of length 3.']);
	end
	k=k+1;
else
	B=obj.MagneticField;
end
if k<nargin
	error('FluxQon:Ion:Hamiltonian:WrongNargin',...
		'Incorrect number of input arguments.');
end

H=kron(diag(obj.MultipletEnergy),eye((2*obj.ElectronSpin+1) ...
	*(2*obj.NuclearSpin+1)))+obj.ZeemanHamiltonian('S',B) ...
	+obj.HyperfineHamiltonian('S')+obj.QuadrupoleHamiltonian('S');
[W,E]=eig(H);
Ed=diag(E);
obj.Energy=Ed;
obj.Eigenstate=W;

if b~='S'
	H=E;
end

if b=='F' && obj.NumIons>1
	Ed=Ed-Ed(1);
	Ed(1)=[];
	LD=numel(Ed);
	FD=obj.FockDimension;
	d=FD*ones(1,LD);
	N=Number(FD);
	H=zeros(FD^LD);
	for k=1:LD
		H=H+Ed(k)*Operator.kron(d,k,N);
	end
end

end
