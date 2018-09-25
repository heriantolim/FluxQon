classdef Circle < Object2D
%% Circle
%
% Constructor syntax:
%  obj=Circle() creates a Circle object with default properties.
%
%  obj=Circle(R) creates a Circle object with side radius R.
%
%  obj=Circle(R,[x,y]) additionally specifies the position of the object.
%
%  obj=Circle(R,[x,y],a) additionally specifies the rotation of the object.
%
%  obj=Circle('PropertyName',PropertyValue,...) sets the properties of the
%    object in Name-Value pair syntax.
%
% Requires package:
%  - MatCommon_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: Object2D.
%
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 24/06/2017
% Last modified: 24/06/2017

properties
	Radius=0;
end

properties (Dependent=true, SetAccess=private)
	Area
	Circumference
end

methods
	function obj=Circle(varargin)
		if nargin==0
			% nothing is true
		elseif isstringscalar(varargin{1})
			n=1;
			while n<nargin
				mp=findprop(obj,varargin{n});
				if isempty(mp)
					error('FluxQon:Object2D:Circle:InvalidInput',...
						'''%s'' is not a property of the Circle class.',varargin{n});
				elseif strcmpi(mp.SetAccess,'public')
					obj.(varargin{n})=varargin{n+1};
				else
					error('FluxQon:Object2D:Circle:DeniedAccess',...
						'Denied access to set the ''%s'' property.',varargin{n});
				end
				n=n+2;
			end
			if n<=nargin
				warning('FluxQon:Object2D:Circle:IgnoredInput',...
					'An extra input argument was provided, but is ignored.');
			end
		elseif nargin==1
			obj.Radius=varargin{1};
		elseif nargin==2
			obj.Radius=varargin{1};
			obj.Position=varargin{2};
		elseif nargin==3
			obj.Radius=varargin{1};
			obj.Position=varargin{2};
			obj.Rotation=varargin{3};
		else
			error('FluxQon:Object2D:Circle:WrongNargin',...
				'Incorrect number of input arguments.');
		end
	end

	function set.Radius(obj,x)
		if isrealscalar(x) && x>=0
			obj.Radius=x;
		else
			error('FluxQon:Object2D:Circle:setRadius:InvalidInput',...
				'Input to set the radius must be a positive real scalar.');
		end
	end
	
	function x=get.Area(obj)
		x=pi*obj.Radius^2;
	end
	
	function x=get.Circumference(obj)
		x=2*pi*obj.Radius;
	end
	
	G=GCF(obj,x,y,z)
end

methods (Static=true)
	Gr=GCFr(R,r,z)
	Gz=GCFz(R,r,z)
end

end
