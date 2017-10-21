function H=round2herm(H,varargin)
%% Round to Hermitian
%  H=obj.round2herm(H) rounds the matrix H until it is hermitian. If the
%    rounding exceeds the relative tolerance of the object then an error is
%	  returned.
%
%  H=obj.round2herm(H,tol) sets the maximum relative tolerance to tol.
%
% Outputs:
%  H   : The rounded matrix, returned as a Hermitian matrix.
%
% Inputs:
%  H   : The matrix to be rounded, specified as a square complex matrix.
%  tol : Tolerance, specified as a real scalar between 0 and 1.
%
% Requires package:
%  - Common_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 25/09/2016
% Last modified: 25/09/2016

% Defaults
REL_TOL=1e-8;

assert(iscomplexmatrix(H) && diff(size(H))==0,...
	'QuantMech:Operator:round2herm:InvalidInput',...
	'The first argument must be a square complex matrix');
assert(nargin<3,...
	'QuantMech:Operator:round2herm:WrongNargin',...
	'Incorrect number of input arguments.');
if nargin>1
	if isrealscalar(varargin{1}) && varargin{1}>0 && varargin{1}<1
		REL_TOL=varargin{1};
	else
		error('QuantMech:Operator:round2herm:InvalidInput',...
			'Input to the tolerance must be a positive real scalar less than 1.');
	end
end

if ishermitian(H)
	return
end
n=ceil(log10(norm(H-H')));
N=floor(log10(norm(H)*REL_TOL));
H=round(H,-n);
while ~ishermitian(H) && n<N
	n=n+1;
	H=round(H,-n);
end
if n>=N
	error('QuantMech:Operator:round2herm:ToleranceExceeded',...
		'The matrix non-hermiticity exceeds the tolerance.');
end

end
