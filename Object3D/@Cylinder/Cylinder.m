classdef Cylinder < Object3D
%% Cylinder
%
% Constructor syntax:
%  obj=Cylinder() creates a Cylinder object with default properties.
%
%  obj=Cylinder(R,H) creates a Cylinder object with side radius R and height H.
%
%  obj=Cylinder(R,H,[x,y,z]) additionally specifies the position of the object.
%
%  obj=Cylinder(R,H,[x,y,z],[a,b,c]) additionally specifies the rotation of the
%    object.
%
%  obj=Cylinder('PropertyName',PropertyValue,...) sets the properties of the
%    object in Name-Value pair syntax.
%
% Requires package:
%  - MatCommon_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: Object3D.
%
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 15/06/2017
% Last modified: 15/06/2017

properties
	Radius=0;
	Height=0;
end

properties (Dependent=true, SetAccess=private)
	Volume
end

methods
	function obj=Cylinder(varargin)
		if nargin==0
			% nothing is true
		elseif isstringscalar(varargin{1})
			n=1;
			while n<nargin
				mp=findprop(obj,varargin{n});
				if isempty(mp)
					error('FluxQon:Object3D:Cylinder:InvalidInput',...
						'''%s'' is not a property of the Cylinder class.',...
						varargin{n});
				elseif strcmpi(mp.SetAccess,'public')
					obj.(varargin{n})=varargin{n+1};
				else
					error('FluxQon:Object3D:Cylinder:DeniedAccess',...
						'Denied access to set the ''%s'' property.',varargin{n});
				end
				n=n+2;
			end
			if n<=nargin
				warning('FluxQon:Object3D:Cylinder:IgnoredInput',...
					'An extra input argument was provided, but is ignored.');
			end
		elseif nargin==2
			obj.Radius=varargin{1};
			obj.Height=varargin{2};
		elseif nargin==3
			obj.Radius=varargin{1};
			obj.Height=varargin{2};
			obj.Position=varargin{3};
		elseif nargin==4
			obj.Radius=varargin{1};
			obj.Height=varargin{2};
			obj.Position=varargin{3};
			obj.Rotation=varargin{4};
		else
			error('FluxQon:Object3D:Cylinder:WrongNargin',...
				'Incorrect number of input arguments.');
		end
	end

	function set.Radius(obj,x)
		if isrealscalar(x) && x>=0
			obj.Radius=x;
		else
			error('FluxQon:Object3D:Cylinder:setRadius:InvalidInput',...
				'Input to set the radius must be a positive real scalar.');
		end
	end
	
	function set.Height(obj,x)
		if isrealscalar(x) && x>=0
			obj.Height=x;
		else
			error('FluxQon:Object3D:Cylinder:setHeight:InvalidInput',...
				'Input to set the height must be a positive real scalar.');
		end
	end
	
	function x=get.Volume(obj)
		x=pi*obj.Radius^2*obj.Height;
	end
end

end
