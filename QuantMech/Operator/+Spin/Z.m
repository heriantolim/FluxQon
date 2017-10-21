function Sz=Z(S)
%% Spin Z Operator
%  Note: hbar is set to 1.
%
%  Sz=Spin.Z(S) returns the spin Z operator for a given total spin S in natural
%  units (hbar=1). S must be a positive half-integer scalar.
%
% Requires package:
%  - Common_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: X, Y, P, M.
% 
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 30/11/2015
% Last modified: 20/05/2017

%% Main
assert(isintegerscalar(2*S) && S>=0,...
	'QuantMech:Operator:Spin:Z:InvalidInput',...
	'Input to the total spin must be a positive half-integer scalar.');
Sz=diag(-S:S);

end