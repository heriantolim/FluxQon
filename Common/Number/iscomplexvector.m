function tf=iscomplexvector(x)
%% Is Input a Vector of Complex Numbers?
%
% See also: iscomplexscalar, iscomplexmatrix, iscomplexarray.
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 14/09/2015
% Last modified: 14/09/2015

tf=isnumeric(x) && isvector(x);

end
