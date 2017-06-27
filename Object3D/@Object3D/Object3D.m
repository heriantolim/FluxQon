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
%    y-axis, and finally by an angle c about the z-axis..
%
%  obj=Object3D('PropertyName',PropertyValue,...) sets the properties of the
%    object in Name-Value pair syntax.
%
% Requires package:
%  - Common_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: Object2D, rotate.
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 06/12/2015
% Last modified: 24/06/2017

properties
	CoordSys=1;
	Position=[0;0;0];
	Rotation=eye(3);
end

properties (Constant=true, Access=protected)
	CoordsSysList={'Catersian','Cylindrical','Spherical'};
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
	
	function set.CoordSys(obj,x)
		L=obj.CoordSysList;
		if isintegerscalar(x) && x>0 && x<numel(L)
			obj.CoordSys=x;
		elseif isstringscalar(x)
			x=strcmpi(x,L);
			if any(x)
				obj.CoordSys=find(x,1,'first');
			else
				error('FluxQon:Object3D:setCoordSys:UnrecognizedInput',...
					'The specified coordinate system is not recognized.');
			end
		else
			error('FluxQon:Object3D:setCoordSys:InvalidInput',...
				['Input to set the coordinate system must be either ',...
					'an integer between 1 and %d, or a string scalar.'],numel(L));
		end
	end

	function set.Position(obj,x)
		if isrealvector(x) && numel(x)==3
			obj.Position=reshape(x,3,1);
		else
			error('FluxQon:Object3D:setPosition:InvalidInput',...
				'Input to set the position must be a real vector of length 3.');
		end
		obj.afterSetPosition;
	end

	function set.Rotation(obj,x)
		try
			x=obj.rotate(x);
		catch ME1
			if isempty(regexpi(ME1.identifier,...
					':WrongNargin$|:InvalidInput$','once'))
				rethrow(ME1);
			else
				ME=MException('FluxQon:Object3D:setRotation:InvalidInput',...
					'Input to set the rotation failed validation.');
				ME=addCause(ME,ME1);
				throw(ME);
			end
		end
		obj.Rotation=x;
		obj.afterSetRotation;
	end
	
	R=rotate(obj,varargin)
end

methods (Access=protected)
	function afterSetPosition(obj) %#ok<MANU>
	end

	function afterSetRotation(obj) %#ok<MANU>
	end
end

end
