function Sx=X(S)
%% Spin X Operator
%  Note: hbar is set to 1.
%
%  Sz=Spin.X(S) returns the spin X operator for a given total spin S in natural
%  units (hbar=1). S must be a positive half-integer scalar.
%
% Requires package:
%  - Common_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: Y, Z, P, M.
% 
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 30/11/2015
% Last modified: 30/11/2015

%% Main
Sx=(Spin.P(S)+Spin.M(S))/2;

end