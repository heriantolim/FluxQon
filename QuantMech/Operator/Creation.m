function A=Creation(N)
%% Creation Operator
%  A=Creation(N) returns the matrix representation for the creation operator.
%  The matrix is truncated to the N lowest Fock states.
%
% See also: Annihilation, Number.
% 
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 20/05/2017
% Last modified: 20/05/2017

A=diag(sqrt(1:N-1),-1);

end