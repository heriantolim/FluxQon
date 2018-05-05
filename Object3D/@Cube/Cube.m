classdef Cube < Object3D
%% Cube
%
% Constructor syntax:
%  obj=Cube() creates a Cube object with default properties.
%
%  obj=Cube(L) creates a Cube object with side length L.
%
%  obj=Cube(L,[x,y,z]) additionally specifies the position of the object.
%
%  obj=Cube(L,[x,y,z],[a,b,c]) additionally specifies the rotation of the
%    object.
%
%  obj=Cube('PropertyName',PropertyValue,...) sets the properties of the object
%    in Name-Value pair syntax.
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
	Length=0;
end

properties (Dependent=true, SetAccess=private)
	Volume
end

methods
	function obj=Cube(varargin)
		if nargin==0
			% nothing is true
		elseif isstringscalar(varargin{1})
			n=1;
			while n<nargin
				mp=findprop(obj,varargin{n});
				if isempty(mp)
					error('FluxQon:Object3D:Cube:InvalidInput',...
						'''%s'' is not a property of the Cube class.',varargin{n});
				elseif strcmpi(mp.SetAccess,'public')
					obj.(varargin{n})=varargin{n+1};
				else
					error('FluxQon:Object3D:Cube:DeniedAccess',...
						'Denied access to set the ''%s'' property.',varargin{n});
				end
				n=n+2;
			end
			if n<=nargin
				warning('FluxQon:Object3D:Cube:IgnoredInput',...
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
			error('FluxQon:Object3D:Cube:WrongNargin',...
				'Incorrect number of input arguments.');
		end
	end

	function set.Length(obj,x)
		if isrealscalar(x) && x>=0
			obj.Length=x;
		else
			error('FluxQon:Object3D:Cube:setLength:InvalidInput',...
				'Input to set the length must be a positive real scalar.');
		end
	end
	
	function x=get.Volume(obj)
		x=obj.Length^3;
	end
end

end
