function S=vec2trlsymm(v)
%% Vector to Traceless Symmetric Matrix
%  S=vec2trlsymm(v) converts the vector v into a traceless symmetric matrix S.
%  Each element in v is put into the upper triangular block in S in sequence
%  from left to right, down to bottom. The lower triangular block is obtained by
%  taking a transpose of the upper triangular block.
%
% Requires package:
%  - Common_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: vec2herm, vec2symm, vec2trlherm.
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 26/06/2017
% Last modified: 26/06/2017

if iscomplexvector(v)
	n=numel(v);
	n=roots([1,-1,-2*n]);
	n=max(n);
	if floor(n)>=n
		S=zeros(n);
		k=0;
		for i=1:n-1
			for j=i+1:n
				k=k+1;
				S(i,j)=v(k);
				S(j,i)=v(k);
			end
		end
	else
		error('QuantMech:Operator:vec2trlsymm:InvalidInput',...
			'Incorrect number of elements in the input vector.');
	end
else
	error('QuantMech:Operator:vec2trlsymm:InvalidInput',...
		'The input argument must be a complex vector.');
end

end