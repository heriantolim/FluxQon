function [E,W,dX,dY]=solveTISE(obj,varargin)
%% Solve Time Independent Schrodinger Equation
%  This method solves the time independent Schrodinger equation for the
%  eigen-energies and eigen-states using the finite difference method with an
%  eight-order accuracy. The solution is assumed periodic with a period 2*pi.
%
%  The object's energy and tunneling energy will be automatically set after a
%  successful run of this method.
%
% Optional input:
%  - N : Number of points, specified as an integer scalar > 9, defaults to 100.
%
% Outputs:
%  - E  : Eigen-energies, returned as a real row vector of length N^2.
%  - W  : Eigen-states, returned as a real or complex matrix of size NxNxN^2.
%  - dX : Infinitesimal length X.
%  - dY : Infinitesimal length Y.
%
% Requires package:
%  - MatCommon_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 13/07/2016
% Last modified: 15/07/2016

% Retrieve required quantities
EJ=obj.JosephsonEnergy;
ECp=obj.CoulombPlus;
ECm=obj.CoulombMinus;
assert(~isempty(EJ) && ~isempty(ECp) && ~isempty(ECm),...
	'FluxQon:Sq3JJ:solveTISE:MissingData',...
	['Please set the Coulomb and Josephson energies, ',...
		'prior to calling this method.']);

% inputs validation
if nargin==1
	N=100;
elseif nargin==2
	N=varargin{1};
	assert(isintegerscalar(N) && N>9,...
		'FluxQon:Sq3JJ:solveTISE:InvalidInput',...
		'Input to the number of points must be a large-enough integer scalar.');
else
	error('FluxQon:Sq3JJ:solveTISE:WrongNargin',...
		'Incorrect number of input arguments.');
end
PB=obj.PhaseBias;% default to pi
PT=obj.PhaseTune;% default to 0
a=abs(cos(PT/2))*obj.AlphaCoefficient;
N2=N^2;

% scale the orders to avoid errors due to machine epsilon
ECp=ECp/EJ;
ECm=ECm/EJ;

% the finite difference coefficients
FDC=[-205/72,8/5,-1/5,8/315,-1/560];

% matrix construction
X=linspace(-pi/2,pi/2,N).';
dX=X(2)-X(1);
dY=dX;
U=zeros(N,N);
for i=1:N
	for j=1:N
		U(i,j)=-2*cos(X(i))*cos(X(j))-a*cos(2*X(i)+PB);
	end
end
U=reshape(U,1,N2);

F=numel(FDC)-1;
Fp=-ECp/2/dX^2*FDC;
Kp=Fp(1)*eye(N);
for f=1:F
	Kp=Kp+Fp(f+1)*(diag(ones(1,f),f-N)+diag(ones(1,N-f),f) ...
		+diag(ones(1,f),N-f)+diag(ones(1,N-f),-f));
end
Fm=-ECm/2/dY^2*FDC;
Km=Fm(1)*eye(N);
for f=1:F
	Km=Km+Fm(f+1)*(diag(ones(1,f),f-N)+diag(ones(1,N-f),f) ...
		+diag(ones(1,f),N-f)+diag(ones(1,N-f),-f));
end

% eigen-decomposition
[W,E]=eig(kron(Km,Kp)+diag(U));
E=EJ*reshape(diag(E),1,N2);
W=reshape(W,N,N,N2);

% Note that W(x,y) is a normalized set of discrete points, not a continuous
% function it ought to be mathematically. The scale of W(x,y) is dependent on
% the number of points, N. However, the inner-product integration of W(x,y) is a
% constant and should be computed as a sum, without multiplying it by the
% infinitesimal distances, delta x delta y.

% update object with the new results
obj.Energy=E(2)-E(1);
obj.TunnelingEnergy=2*a*EJ ...
	*trapz(trapz(conj(W(:,:,1)).*repmat(sin(2*X+PB),1,N).*W(:,:,2)));

% Renormalize W(x,y) so that its values are independent of N. Now inner-product
% integrations have to include a multiplication by delta x delta y.
W=W/dX;

end