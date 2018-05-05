function [H,d]=Hamiltonian(varargin)
%% Construct Hamiltonian
%  H=Construct.Hamiltonian(obj1,obj2,obj3,...,objN) returns the total
%    Hamiltonian of obj1,obj2,...,objN, including the interactions between the
%    objects.
%
%  H=Construct.Hamiltonian(d,n1,obj1,n2,obj2,...,nN,objN) returns the total
%    Hamiltonian with the kronecker product dimensions specified in d. The
%    Hilbert subspace for obj_i is positioned at index n_i.
%
%  H=Construct.Hamiltonian(obj1,obj2,obj3,...,objN,[])
%  H=Construct.Hamiltonian(d,n1,obj1,n2,obj2,...,nN,objN,[])
%    returns the Hamiltonian without the interactions.
%
%  H=Construct.Hamiltonian(obj1,obj2,obj3,...,objN,CIX)
%  H=Construct.Hamiltonian(d,n1,obj1,n2,obj2,...,nN,objN,CIX)
%    returns the Hamiltonian with the interactions between any two objects in
%    the combination pairs specified in the matrix CIX. CIX is a matrix of size
%    M x 2 where the elements are integers {1,2,3,...,N} that correspond to the
%    indices of the objects in the list. Each row in CIX specifies a combination
%    between two indices (Xi,Xj), which is used to mark that the interaction
%    between obj_i and obj_j should be added into the Hamiltonian.
%
%  H=Construct.Hamiltonian(...,'RWA') uses the RWA to construct the interaction
%    Hamiltonian.
%
% Optional outputs:
%  - d : List of the Hilbert-subspace dimensions, returned as an integer vector.
%        The product of the elements in d is equal to the dimension of H.
%
% Requires package:
%  - MatCommon_v1.0.0+
%
% Tested on:
%  - MATLAB R2017a
%
% See also: Interaction, Lindblad.
%
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 16/06/2017
% Last modified: 01/08/2017

H=0;
d=1;
K=nargin;
if K==0
	return
end

if isstringscalar(varargin{K})
	K=K-1;
end

if isintegermatrix(varargin{K})
	K=K-1;
end

if isintegervector(varargin{1})
	J=(K-1)/2;
	assert(floor(J)>=J,...
		'FluxQon:Construct:Hamiltonian:WrongNargin',...
		'Incorrect number of input arguments.');
	assert(all(varargin{1}>0),...
		'FluxQon:Construct:Hamiltonian:InvalidInput',...
		'Input to the subspace dimensions must be a positive integer vector.');
	d=varargin{1};
	k=2*(1:J);
	assert(all(cellfun(@(x)isintegerscalar(x) && x>0,varargin(k))),...
		'FluxQon:Construct:Hamiltonian:InvalidInput',...
		'Input to the subspace indices must be a positive integer scalar.');
	n=[varargin{k}];
	k=k+1;
else
	J=K;
	d=[];
	n=1:J;
	k=n;
end

assert(all(cellfun(@(x)ismethod(x,'Hamiltonian'),varargin(k))),...
	'FluxQon:Construct:Hamiltonian:InvalidInput',...
	'The input objects must have a method or a property named ''Hamiltonian''.');

if isempty(d)
	d=Hilbert.dimension(varargin{k});
end

for j=1:J
	H=H+Operator.kron(d,n(j),varargin{k(j)}.Hamiltonian);
end

H=H+Construct.Interaction(varargin{:});

end
