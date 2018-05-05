function Sp=P(S,N,F)
%% Holstein-Primakoff Spin Raising Operator
%  Note: hbar is set to 1.
%
%  Sp=Spin.P(S,N,F) returns the approximated, collective spin raising operator
%  of N identical atoms, each with a total spin S. The collective operator is
%  formed from a combination of the annihilation-and-creation operators, each
%  with Fock dimension F, by using the Holstein-Primakoff transformation [1,2].
%
% Input arguments:
%  S : Total spin of each atom, specified as a positive half-integer scalar.
%  N : Number of atoms, specified as an integer scalar greater than (F-1)^(2S).
%  F : Fock dimension of each annihilation and creation operator.
%
% References:
%  [1] T. Holstein and H. Primakoff, Phys. Rev. 58, 1098 (1940).
%  [2] Z. Kurucz and K. Molmer, Phys. Rev. A 81, 1 (2010).
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
% First created: 08/06/2017
% Last modified: 08/06/2017

M=2*S;
assert(isintegerscalar(M) && M>=0,...
	'FluxQon:Operator:HPSpin:Z:InvalidInput',...
	'Input to the total spin must be a positive half-integer scalar.');
assert(isintegerscalar(F) && F>0,...
	'FluxQon:Operator:HPSpin:Z:InvalidInput',...
	'Input to the Fock dimension must be a positive integer scalar.');
assert(isintegerscalar(N) && N>=(F-1)^M,...
	'FluxQon:Operator:HPSpin:Z:InvalidInput',...
	['Input to the number of atoms must be an integer scalar greater ',...
		'or equal to %d.'],(F-1)^M);

if M==0
	Sp=0;
else
	D=F*ones(1,M);
	Sp=N*eye(F^M);
	A=Number(F);
	for m=1:M
		Sp=Sp-Operator.kron(D,m,A);
	end
	A=Annihilation(F);
	C=Creation(F);
	Sp=Operator.kron(D,1,C)*sqrt(M*Sp);
	M=M-1;
	for m=2:M
		Sp=Sp+sqrt(S*(S+1)-(m-S)*(m-S+1)) ...
			*Operator.kron(D,m+1,C,m,A);
	end
end

end
