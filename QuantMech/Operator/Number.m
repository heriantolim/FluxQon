function A=Number(N)
%% Number Operator
%  A=Number(N) returns the matrix representation for the number operator. The
%  matrix is truncated to the N lowest Fock states.
%
% See also: Annihilation, Creation.
% 
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 20/05/2017
% Last modified: 20/05/2017

A=diag(0:N-1);

end