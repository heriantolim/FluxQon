function [E,W,dX]=solveTISE(obj,varargin)
%% Solve Time Independent Schrodinger Equation
%  This method solves the time independent Schrodinger equation for the
%  eigen-energies and eigen-states using the finite difference method with an
%  eight-order accuracy. The solution is assumed periodic with a period 4*pi.
%
%  The object's energy and tunneling energy will be automatically set after a
%  successful run of this method.
%
% Optional input:
%  - N : Number of points, specified as an integer scalar > 9, defaults to 1000.
%
% Outputs:
%  - E  : Eigen-energies, returned as a real row vector of length N.
%  - W  : Eigen-states, returned as a real or complex matrix of size N x N.
%  - dX : Infinitesimal length X.
%
% Requires package:
%  - MatCommon_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 13/12/2015
% Last modified: 15/07/2016

% Retrieve required quantities
EC=obj.CoulombEnergy;
EL=obj.InductanceEnergy;
EJ=obj.JosephsonEnergy;
assert(~(isempty(EC) || isempty(EL) || isempty(EJ)),...
	'FluxQon:Qubit:SQUID:solveTISE:MissingData',...
	['Please set the Coulomb, inductance, and Josephson energies, ',...
		'prior to calling this method.']);

% Inputs validation
if nargin==1
	N=1000;
elseif nargin==2
	N=varargin{1};
	assert(isintegerscalar(N) && N>9,...
		'FluxQon:Qubit:SQUID:solveTISE:InvalidInput',...
		'Input to the number of points must be a large-enough integer scalar.');
else
	error('FluxQon:Qubit:SQUID:solveTISE:WrongNargin',...
		'Incorrect number of input arguments.');
end
PB=obj.PhaseBias;% default to pi
PT=obj.PhaseTune;% default to 0

% scale the orders to avoid errors due to machine epsilon
EC=EC/EL;
EJ=abs(cos(PT/2))*EJ/EL;

% the finite difference coefficients
FDC=[-205/72,8/5,-1/5,8/315,-1/560];

% matrix construction
X=linspace(0,2*pi,N).';
dX=X(2)-X(1);
F=numel(FDC)-1;
H=zeros(N,N);
FDC=-EC/2/dX^2*FDC;
for f=1:F
	H=H+FDC(f+1)*(diag(ones(1,f),f-N)+diag(ones(1,N-f),f)...
		+diag(ones(1,f),N-f)+diag(ones(1,N-f),-f));
end

% The domain X is shifted by pi-PB for better accuracy in calculations where PB
% is close to the boundary.
H=H+diag(FDC(1)+(X-pi).^2/2-EJ*cos(X-(pi-PB)));

% eigen-decomposition
[W,E]=eig(H);
E=EL*reshape(diag(E),1,N);

% Note that W(x) is a normalized set of discrete points, not a continuous
% function it ought to be mathematically. The scale of W(x) is dependent on the
% number of points, N. However, the inner-product integration of W(x) is a
% constant and should be computed as a sum, without multiplying it by the
% infinitesimal distance, delta x.

% update object with the new results
obj.Energy=E(2)-E(1);
obj.TunnelingEnergy=-2*EL*trapz(conj(W(:,1)).*(X-pi).*W(:,2));

% Renormalize W(x) so that its values are independent of N. Now inner-product
% integrations have to include a multiplication by delta x.
W=W/sqrt(dX);

end
