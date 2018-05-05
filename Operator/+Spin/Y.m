function Sy=Y(S)
%% Spin Y Operator
%  Note: hbar is set to 1.
%
%  Sz=Spin.Y(S) returns the spin Y operator for a given total spin S in natural
%  units (hbar=1). S must be a positive half-integer scalar.
%
% Requires package:
%  - MatCommon_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: X, Z, P, M.
% 
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 30/11/2015
% Last modified: 30/11/2015

%% Main
Sy=(Spin.P(S)-Spin.M(S))/2i;

end