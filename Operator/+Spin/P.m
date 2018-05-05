function Sp=P(S)
%% Spin Raising Operator
%  Note: hbar is set to 1.
%
%  Sp=Spin.plus(S) returns the spin raising operator for a given total spin S
%  in natural units (hbar=1). S must be a positive half-integer scalar.
%
% Requires package:
%  - MatCommon_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: X, Y, Z, M.
% 
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 30/11/2015
% Last modified: 08/06/2017

%% Main
S=2*S;
assert(isintegerscalar(S) && S>=0,...
	'FluxQon:Operator:Spin:plus:InvalidInput',...
	'Input to the total spin must be a positive half-integer scalar.');
if S==0
	Sp=0;
else
	Sc=S:-2:-S;
	Sp=zeros(1,S);
	Sp(1)=Sc(1);
	s=1;
	while s<S
		s=s+1;
		Sp(s)=Sp(s-1)+Sc(s);
	end
	Sp=diag(sqrt(Sp),-1);
end

end