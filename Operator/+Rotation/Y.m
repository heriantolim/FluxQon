function R=Y(a)
%% Y Rotation Matrix
%  R=Rotation.Y(a) returns the matrix for a rotation by an angle a (in radians)
%  about the y axis.
%
% Tested on:
%  - MATLAB R2017a
%
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 18/03/2018
% Last modified: 18/03/2018

assert(isrealscalar(a),...
	'FluxQon:Rotation:X:InvalidInput',...
	'Input to the rotation angle must be a real scalar.');

c=cos(a);
s=sin(a);
R=[c,0,s;0,1,0;-s,0,c];

end