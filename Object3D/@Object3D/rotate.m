function R=rotate(obj,varargin)
%% Rotation Matrix
%  R=Object3D.rotate([a,b,c]) or R=Object3D.rotate(a,b,c) returns the matrix for
%  a rotation: first by an angle a about the z-axis, then by an angle b about
%  the y-axis, and finally by an angle c about the z-axis. The angles are in
%  radians.
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

if nargin==2
	assert(isrealvector(varargin{1}) && numel(varargin{1})==3,...
		'Object3D:rotate:InvalidInput',...
		['Input to the vector of rotation angles must be ',...
			'a real vector of length 3.']);
	a=varargin{1};
elseif nargin==4
	assert(all(cellfun(@isrealscalar,varargin)),...
		'Object3D:rotate:InvalidInput',...
		'Input to the rotation angles must be a real scalar.');
	a=[varargin{:}];
else
	error('Object3D:rotate:WrongNargin',...
		'Incorrect number of input arguments.');
end

switch obj.CoordSys
	case 1% Catersian
		c=cos(a);
		s=sin(a);
		R=[c(3),-s(3),0;s(3),c(3),0;0,0,1] ...
			*[c(2),0,s(2);0,1,0;-s(2),0,c(2)] ...
			*[c(1),-s(1),0;s(1),c(1),0;0,0,1];
	case 2% Cylindrical
		error('Object3D:rotate:IncompleteCode',...
			'Incomplete code.');
	case 3% Spherical
		error('Object3D:rotate:IncompleteCode',...
			'Incomplete code.');
end

end