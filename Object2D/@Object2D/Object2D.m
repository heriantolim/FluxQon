classdef Object2D < handle
%% Object in 2D Space
%
% Constructor syntax:
%  obj=Object2D() creates an Object2D with default properties.
%
%  obj=Object2D([x,y]) creates an Object2D with Position set to [x,y].
%
%  obj=Object2D([x,y],a) additionally sets the Rotation of the object to
%    a rotation by an angle a.
%
%  obj=Object2D('PropertyName',PropertyValue,...) sets the properties of the
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
% See also: Object3D, rotate.
%
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 15/06/2017
% Last modified: 02/04/2018

properties
	CoordinateSystem=1;
	Position=[0;0];
	Rotation=0;
end

properties (Dependent=true)
	RotationMatrix
end

properties (Constant=true, Access=protected)
	CoordinateSystemList={'Catersian','Polar'};
end

methods
	function obj=Object2D(varargin)
		if nargin==0
			% nothing is true
		elseif isstringscalar(varargin{1})
			n=1;
			while n<nargin
				mp=findprop(obj,varargin{n});
				if isempty(mp)
					error('FluxQon:Object2D:InvalidInput',...
						'''%s'' is not a property of the Object2D class.',...
						varargin{n});
				elseif strcmpi(mp.SetAccess,'public')
					obj.(varargin{n})=varargin{n+1};
				else
					error('FluxQon:Object2D:DeniedAccess',...
						'Denied access to set the ''%s'' property.',varargin{n});
				end
				n=n+2;
			end
			if n<=nargin
				warning('FluxQon:Object2D:IgnoredInput',...
					'An extra input argument was provided, but is ignored.');
			end
		elseif nargin==1
			obj.Position=varargin{1};
		elseif nargin==2
			obj.Position=varargin{1};
			obj.Rotation=varargin{2};
		else
			error('FluxQon:Object2D:WrongNargin',...
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
				error('FluxQon:Object2D:setCoordinateSystem:UnrecognizedInput',...
					'The specified coordinate system is not recognized.');
			end
		else
			error('FluxQon:Object2D:setCoordinateSystem:InvalidInput',...
				['Input to set the coordinate system must be either ',...
					'an integer between 1 and %d, or a string scalar.'],numel(L));
		end
	end

	function set.Position(obj,x)
		if isrealvector(x) && numel(x)==2
			obj.Position=x(:);
			obj.afterSetPosition;
		else
			error('FluxQon:Object2D:setPosition:InvalidInput',...
				'Input to set the position must be a real vector of length 2.');
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
		if a==0
			R=eye(2);
		else
			switch obj.CoordinateSystem
				case 1% Catersian
					c=cos(a);
					s=sin(a);
					R=[c,-s;s,c];
				case 2% Polar
					error('FluxQon:Object2D:getRotationMatrix:IncompleteCode',...
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
