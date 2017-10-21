function tf=isintegermatrix(x)
%% Is Input a Matrix of Integer Numbers?
%
% See also: isintegerscalar, isintegervector, isintegerarray.
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 01/03/2013
% Last modified: 01/03/2013

tf=isintegerarray(x) && ismatrix(x);

end
