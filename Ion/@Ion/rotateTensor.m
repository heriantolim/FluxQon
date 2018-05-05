function T=rotateTensor(obj,T)
%% Rotate a Tensor of Order 2 (a 3x3 Matrix)
%  T=obj.rotateTensor(T) rotate a 3x3 matrix according to the rotation matrix
%  defined in the object properties.
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 05/12/2015
% Last modified: 25/09/2016

if ~isscalar(T)
	R=obj.RotationMatrix;
	T=R*T*R.';
end

end
