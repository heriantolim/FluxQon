function H=vec2herm(v)
%% Vector to Hermitian Matrix
%  H=vec2herm(v) converts the vector v into a hermitian matrix H. Each element
%  in v is put into the upper off-diagonal block in H in sequence from left to
%  right, down to bottom. The lower off-diagonal block is obtained by taking a
%  complex conjugate transpose of the upper off-diagonal block.
%
% Requires package:
%  - MatCommon_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: vec2symm, vec2trlherm, vec2trlsymm.
%
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 26/06/2017
% Last modified: 26/06/2017

if iscomplexvector(v)
	n=max(roots([1,1,-2*numel(v)]));
	if floor(n)>=n
		H=zeros(n);
		k=0;
		for i=1:n
			for j=i:n
				k=k+1;
				H(i,j)=v(k);
				H(j,i)=v(k)';
			end
		end
	else
		error('FluxQon:Operator:vec2herm:InvalidInput',...
			'Incorrect number of elements in the input vector.');
	end
else
	error('FluxQon:Operator:vec2herm:InvalidInput',...
		'The input argument must be a complex vector.');
end

end