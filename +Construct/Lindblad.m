function [L,d]=Lindblad(varargin)
%% Construct Lindblad Superoperator
%  L=Construct.Lindblad(obj1,obj2,obj3,...,objN) returns the Lindblad
%    superoperator of obj1,obj2,...,objN, which takes into account the
%    interactions between the objects.
%
%  L=Construct.Lindblad(d,n1,obj1,n2,obj2,...,nN,objN) returns the Lindblad
%    superoperator with the kronecker product dimensions specified in d. The
%    Hilbert subspace for obj_i is positioned at index n_i.
%
%  L=Construct.Lindblad(obj1,obj2,obj3,...,objN,[])
%  L=Construct.Lindblad(d,n1,obj1,n2,obj2,...,nN,objN,[])
%    returns the Lindblad superoperator exclusive of the interaction parts.
%
%  L=Construct.Lindblad(obj1,obj2,obj3,...,objN,CIX)
%  L=Construct.Lindblad(d,n1,obj1,n2,obj2,...,nN,objN,CIX)
%    returns the Lindblad superoperator with the interactions between any two
%    objects in the combination pairs specified in the matrix CIX. CIX is a
%    matrix of size M x 2 where the elements are integers {1,2,3,...,N} that
%    correspond to the indices of the objects in the list. Each row in CIX
%    specifies a combination between two indices (Xi,Xj), which is used to mark
%    that the interaction between obj_i and obj_j should be added into the
%    calculation of L.
%
% Optional outputs:
%  - d : List of the Hilbert-subspace dimensions, returned as an integer vector.
%        The product of the elements in d is equal to the square-root dimension
%        of L.
%
% Requires package:
%  - Common_v1.0.0+
%  - QuantMech_v1.0.0+
%
% Tested on:
%  - MATLAB R2017a
%
% See also: Hamiltonian, Interaction.
% 
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 19/06/2017
% Last modified: 19/06/2017

K=nargin;
if K==0
	L=0;
	d=1;
	return
end

if isintegermatrix(varargin{K})
	K=K-1;
end

if isintegervector(varargin{1})
	J=(K-1)/2;
	assert(floor(J)>=J,...
		'FluxQon:Construct:Lindblad:WrongNargin',...
		'Incorrect number of input arguments.');
	assert(all(varargin{1}>0),...
		'FluxQon:Construct:Lindblad:InvalidInput',...
		'Input to the subspace dimensions must be a positive integer vector.');
	d=varargin{1};
	k=2*(1:J);
	assert(all(cellfun(@(x)isintegerscalar(x) && x>0,varargin(k))),...
		'FluxQon:Construct:Lindblad:InvalidInput',...
		'Input to the subspace indices must be a positive integer scalar.');
	n=[varargin{k}];
	k=k+1;
else
	J=K;
	d=[];
	n=1:J;
	k=n;
end

assert(all(cellfun(@(x)ismethod(x,'Lindblad'),varargin(k))),...
	'FluxQon:Construct:Lindblad:InvalidInput',...
	'The input objects must have a method or a property named ''Lindblad''.');

if isempty(d)
	d=Hilbert.dimension(varargin{k});
end

L=Construct.Interaction(varargin{:});
if ~isscalar(L)
	D=prod(d);
	L=L/Constant.ReducedPlanck;
	L=1i*(kron(L',eye(D))-kron(eye(D),L));
end

for j=1:J
	L=L+varargin{k(j)}.Lindblad(d,n(j));
end

end