classdef Square < Object2D
%% Square
%
% Constructor syntax:
%  obj=Square() creates a Square object with default properties.
%
%  obj=Square(L) creates a Square object with side length L.
%
%  obj=Square(L,[x,y]) additionally specifies the position of the object.
%
%  obj=Square(L,[x,y],a) additionally specifies the rotation of the object.
%
%  obj=Square('PropertyName',PropertyValue,...) sets the properties of the
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
% First created: 15/06/2017
% Last modified: 15/06/2017

properties
	Length=0;
end

properties (Dependent=true, SetAccess=private)
	Area
	Circumference
end

methods
	function obj=Square(varargin)
		if nargin==0
			% nothing is true
		elseif isstringscalar(varargin{1})
			n=1;
			while n<nargin
				mp=findprop(obj,varargin{n});
				if isempty(mp)
					error('FluxQon:Object2D:Square:InvalidInput',...
						'''%s'' is not a property of the Square class.',varargin{n});
				elseif strcmpi(mp.SetAccess,'public')
					obj.(varargin{n})=varargin{n+1};
				else
					error('FluxQon:Object2D:Square:DeniedAccess',...
						'Denied access to set the ''%s'' property.',varargin{n});
				end
				n=n+2;
			end
			if n<=nargin
				warning('FluxQon:Object2D:Square:IgnoredInput',...
					'An extra input argument was provided, but is ignored.');
			end
		elseif nargin==1
			obj.Length=varargin{1};
		elseif nargin==2
			obj.Length=varargin{1};
			obj.Position=varargin{2};
		elseif nargin==3
			obj.Length=varargin{1};
			obj.Position=varargin{2};
			obj.Rotation=varargin{3};
		else
			error('FluxQon:Object2D:Square:WrongNargin',...
				'Incorrect number of input arguments.');
		end
	end

	function set.Length(obj,x)
		if isrealscalar(x) && x>=0
			obj.Length=x;
		else
			error('FluxQon:Object2D:Square:setLength:InvalidInput',...
				'Input to set the length must be a positive real scalar.');
		end
	end
	
	function x=get.Area(obj)
		x=obj.Length^2;
	end
	
	function x=get.Circumference(obj)
		x=4*obj.Length;
	end
end

end
