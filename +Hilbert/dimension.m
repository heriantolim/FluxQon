function d=dimension(varargin)
%% Hilbert Dimension
%  d=Hilbert.dimension(obj1,obj2,...) returns the Hilbert dimensions of each
%  objects in obj1 array, obj2 array, and so on.
%
% Requires package:
%  - MatCommon_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
% 
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 18/12/2015
% Last modified: 18/12/2015

d=[];
N=nargin;
for n=1:N
	assert(all(isprop(varargin{n},'HilbertDimension')),...
		'FluxQon:Hilbert:dimension:InvalidInput',...
		'The input objects must have the propery: HilbertDimension.');
	obj=varargin{n}(:);
	O=numel(obj);
	for o=1:O
		d=[d,obj(o).HilbertDimension]; %#ok<AGROW>
	end
end

end