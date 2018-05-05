function Sm=M(S)
%% Spin Lowering Operator
%  Note: hbar is set to 1.
%
%  Sp=Spin.minus(S) returns the spin lowering operator for a given total spin S
%  in natural units (hbar=1). S must be a positive half-integer scalar.
%
% Requires package:
%  - MatCommon_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: X, Y, Z, P.
% 
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 30/11/2015
% Last modified: 08/06/2017

%% Main
S=2*S;
assert(isintegerscalar(S) && S>=0,...
	'FluxQon:Operator:Spin:minus:InvalidInput',...
	'Input to the total spin must be a positive half-integer scalar.');
if S==0
	Sm=0;
else
	Sc=S:-2:-S;
	Sm=zeros(1,S);
	Sm(1)=Sc(1);
	s=1;
	while s<S
		s=s+1;
		Sm(s)=Sm(s-1)+Sc(s);
	end
	Sm=diag(sqrt(Sm),1);
end

end