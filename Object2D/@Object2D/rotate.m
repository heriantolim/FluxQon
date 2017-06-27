function R=rotate(obj,a)
%% Rotation Matrix
%  R=Object3D.rotate(a) returns the matrix for a rotation by an angle a (in
%  radians).
%
% Output:
%  R: a 3x3 rotation matrix.
%
% Requires package:
%  - Common_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 24/09/2016
% Last modified: 24/06/2017

assert(isrealscalar(a),...
	'FluxQon:Object2D:rotate:InvalidInput',...
	'Input to the rotation angle must be a real scalar.');

switch obj.CoordSys
	case 1% Catersian
		c=cos(a);
		s=sin(a);
		R=[c,-s;s,c];
	case 2% Polar
		error('Object2D:rotate:IncompleteCode',...
			'Incomplete code.');
end

end