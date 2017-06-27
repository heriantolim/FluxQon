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
%  - Common_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: Object3D, rotate.
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 15/06/2017
% Last modified: 24/06/2017

properties
	CoordSys=1;
	Position=[0;0];
	Rotation=eye(2);
end

properties (Constant=true, Access=protected)
	CoordsSysList={'Catersian','Polar'};
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
	
	function set.CoordSys(obj,x)
		L=obj.CoordSysList;
		if isintegerscalar(x) && x>0 && x<numel(L)
			obj.CoordSys=x;
		elseif isstringscalar(x)
			x=strcmpi(x,L);
			if any(x)
				obj.CoordSys=find(x,1,'first');
			else
				error('FluxQon:Object2D:setCoordSys:UnrecognizedInput',...
					'The specified coordinate system is not recognized.');
			end
		else
			error('FluxQon:Object2D:setCoordSys:InvalidInput',...
				['Input to set the coordinate system must be either ',...
					'an integer between 1 and %d, or a string scalar.'],numel(L));
		end
	end

	function set.Position(obj,x)
		if isrealvector(x) && numel(x)==2
			obj.Position=reshape(x,2,1);
		else
			error('FluxQon:Object2D:setPosition:InvalidInput',...
				'Input to set the position must be a real vector of length 2.');
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
				ME=MException('FluxQon:Object2D:setRotation:InvalidInput',...
					'Input to set the rotation failed validation.');
				ME=addCause(ME,ME1);
				throw(ME);
			end
		end
		obj.Rotation=x;
		obj.afterSetRotation;
	end
	
	R=rotate(obj,a)
end

methods (Access=protected)
	function afterSetPosition(obj) %#ok<MANU>
	end

	function afterSetRotation(obj) %#ok<MANU>
	end
end

end
