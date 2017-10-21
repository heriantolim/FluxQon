function Sx=X(S,N,F)
%% Holstein-Primakoff Spin X Operator
%  Note: hbar is set to 1.
%
%  Sx=Spin.X(S,N,F) returns the approximated, collective spin X operator of N
%  identical atoms, each with a total spin S. The collective operator is
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
% First created: 08/06/2017
% Last modified: 08/06/2017

%% Main
Sx=(HPSpin.P(S,N,F)+HPSpin.M(S,N,F))/2;

end
