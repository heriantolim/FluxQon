function H=vec2trlherm(v)
%% Vector to Traceless Hermitian Matrix
%  H=vec2trlherm(v) converts the vector v into a traceless hermitian matrix H.
%  Each element in v is put into the upper off-diagonal block in H in sequence
%  from left to right, down to bottom. The lower off-diagonal block is obtained
%  by taking a complex conjugate transpose of the upper off-diagonal block.
%
% Requires package:
%  - MatCommon_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: vec2herm, vec2symm, vec2trlsymm.
%
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 26/06/2017
% Last modified: 26/06/2017

if iscomplexvector(v)
	n=numel(v);
	n=roots([1,-1,-2*n]);
	n=max(n);
	if floor(n)>=n
		H=zeros(n);
		k=0;
		for i=1:n-1
			for j=i+1:n
				k=k+1;
				H(i,j)=v(k);
				H(j,i)=v(k)';
			end
		end
	else
		error('FluxQon:Operator:vec2trlherm:InvalidInput',...
			'Incorrect number of elements in the input vector.');
	end
else
	error('FluxQon:Operator:vec2trlherm:InvalidInput',...
		'The input argument must be a complex vector.');
end

end