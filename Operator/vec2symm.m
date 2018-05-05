function S=vec2symm(v)
%% Vector to Symmetric Matrix
%  S=vec2symm(v) converts the vector v into a symmetric matrix S. Each element
%  in v is put into the upper triangular block in S in sequence from left to
%  right, down to bottom. The lower triangular block is obtained by taking a
%  transpose of the upper triangular block.
%
% Requires package:
%  - MatCommon_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: vec2herm, vec2trlherm, vec2trlsymm.
%
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 26/06/2017
% Last modified: 26/06/2017

if iscomplexvector(v)
	n=numel(v);
	n=roots([1,1,-2*n]);
	n=max(n);
	if floor(n)>=n
		S=zeros(n);
		k=0;
		for i=1:n
			for j=i:n
				k=k+1;
				S(i,j)=v(k);
				S(j,i)=v(k);
			end
		end
	else
		error('FluxQon:Operator:vec2symm:InvalidInput',...
			'Incorrect number of elements in the input vector.');
	end
else
	error('FluxQon:Operator:vec2symm:InvalidInput',...
		'The input argument must be a complex vector.');
end

end