classdef Object3D < handle
%% Object in 3D Space
%
% Constructor syntax:
%  obj=Object3D() creates an Object3D with default properties.
%
%  obj=Object3D([x,y,z]) creates an Object3D with Position set to [x,y,z].
%
%  obj=Object3D([x,y,z],[a,b,c]) additionally sets the Rotation of the object
%    to a rotation by an angle a about the z-axis, then by an angle b about the
%    y-axis, and finally by an angle c about the z-axis.
%
%  obj=Object3D('PropertyName',PropertyValue,...) sets the properties of the
%    object in Name-Value pair syntax.
%
% Requires package:
%  - MatCommon_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%  - MATLAB R2018a
%
% See also: Object2D, rotate.
%
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 06/12/2015
% Last modified: 02/04/2018

properties
	CoordinateSystem=1;
	Position=[0;0;0];
	Rotation=[0;0;0];
end

properties (Dependent=true)
	RotationMatrix
end

properties (Constant=true, Access=protected)
	CoordinateSystemList={'Catersian','Cylindrical','Spherical'};
end

methods
	function obj=Object3D(varargin)
		if nargin==0
			% nothing is true
		elseif isstringscalar(varargin{1})
			n=1;
			while n<nargin
				mp=findprop(obj,varargin{n});
				if isempty(mp)
					error('FluxQon:Object3D:InvalidInput',...
						'''%s'' is not a property of the Object3D class.',...
						varargin{n});
				elseif strcmpi(mp.SetAccess,'public')
					obj.(varargin{n})=varargin{n+1};
				else
					error('FluxQon:Object3D:DeniedAccess',...
						'Denied access to set the ''%s'' property.',varargin{n});
				end
				n=n+2;
			end
			if n<=nargin
				warning('FluxQon:Object3D:IgnoredInput',...
					'An extra input argument was provided, but is ignored.');
			end
		elseif nargin==1
			obj.Position=varargin{1};
		elseif nargin==2
			obj.Position=varargin{1};
			obj.Rotation=varargin{2};
		else
			error('FluxQon:Object3D:WrongNargin',...
				'Incorrect number of input arguments.');
		end
	end
	
	function set.CoordinateSystem(obj,x)
		L=obj.CoordinateSystemList;
		if isintegerscalar(x) && x>0 && x<=numel(L)
			obj.CoordinateSystem=x;
		elseif isstringscalar(x)
			x=strcmpi(x,L);
			if any(x)
				obj.CoordinateSystem=find(x,1,'first');
			else
				error('FluxQon:Object3D:setCoordinateSystem:UnrecognizedInput',...
					'The specified coordinate system is not recognized.');
			end
		else
			error('FluxQon:Object3D:setCoordinateSystem:InvalidInput',...
				['Input to set the coordinate system must be either ',...
					'an integer between 1 and %d, or a string scalar.'],numel(L));
		end
	end

	function set.Position(obj,x)
		if isrealvector(x) && numel(x)==3
			obj.Position=x(:);
			obj.afterSetPosition;
		else
			error('FluxQon:Object3D:setPosition:InvalidInput',...
				'Input to set the position must be a real vector of length 3.');
		end
	end

	function set.Rotation(obj,x)
		if isrealvector(x) && numel(x)==3
			obj.Rotation=mod(x(:),2*pi);
			obj.afterSetRotation;
		else
			error('FluxQon:Object3D:setRotation:InvalidInput',...
				'Input to set the position must be a real vector of length 3.');
		end
	end
	
	function R=get.RotationMatrix(obj)
		a=obj.Rotation;
		if all(a==0)
			R=eye(3);
		else
			switch obj.CoordinateSystem
				case 1% Catersian
					c=cos(a);
					s=sin(a);
					R=[c(3),-s(3),0;s(3),c(3),0;0,0,1] ...
						*[c(2),0,s(2);0,1,0;-s(2),0,c(2)] ...
						*[c(1),-s(1),0;s(1),c(1),0;0,0,1];
				case 2% Cylindrical
					error('FluxQon:Object3D:getRotationMatrix:IncompleteCode',...
						'Incomplete code.');
				case 3% Spherical
					error('FluxQon:Object3D:getRotationMatrix:IncompleteCode',...
						'Incomplete code.');
			end
		end
	end
end

methods (Access=protected)
	function afterSetPosition(obj) %#ok<MANU>
	end

	function afterSetRotation(obj) %#ok<MANU>
	end
end

end
