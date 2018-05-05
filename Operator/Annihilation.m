function A=Annihilation(N)
%% Annihilation Operator
%  A=Annihilation(N) returns the matrix representation for the annihilation
%  operator. The matrix is truncated to the N lowest Fock states.
%
% See also: Creation, Number.
% 
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 20/05/2017
% Last modified: 20/05/2017

A=diag(sqrt(1:N-1),1);

end