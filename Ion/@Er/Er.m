classdef Er < Ion
%% Er3+ Ion in a Microwave Field
%
% Constructor syntax:
%  obj=Er()
%  obj=Er(0)
%    creates an Er object with a random isotope according to its natural
%    abundance: 22.95% chance being odd, and the rest being even.
%
%  obj=Er(1)
%  obj=Er('odd')
%    creates an Er object with an odd isotope.
%
%  obj=Er(2)
%  obj=Er('even')
%    creates an Er object with an even isotope.
%
%  obj=Er(N,F,...), where N>2, creates an Er object representing a collection of
%    N Er ions. F is the dimension of the annihilation-and-creation operators
%    used to approximate the spin operators.
%
%  obj=Er('PropertyName',PropertyValue,...) sets the properties of the object in
%    Name-Value pair syntax.
%
% Requires package:
%  - MatCommon_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: Ion.
%
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 05/12/2015
% Last modified: 17/06/2017

methods
	function obj=Er(varargin)
		N=nargin;
		if N>1 && isintegerscalar(varargin{1}) && varargin{1}>2
			obj.NumIons=varargin{1};
			obj.FockDimension=varargin{2};
			varargin=varargin(3:N);
			N=N-2;
		end
		if N==0
			obj.Isotope=0;
		elseif isstringscalar(varargin{1})
			n=1;
			while n<N
				mp=findprop(obj,varargin{n});
				if isempty(mp)
					error('FluxQon:Ion:Er:InvalidInput',...
						'''%s'' is not a property of the Er class.',varargin{n});
				elseif strcmpi(mp.SetAccess,'public')
					obj.(varargin{n})=varargin{n+1};
				else
					error('FluxQon:Ion:Er:DeniedAccess',...
						'Denied access to set the ''%s'' property.',varargin{n});
				end
				n=n+2;
			end
			if n<=N
				warning('FluxQon:Ion:Er:IgnoredInput',...
					'An extra input argument was provided, but is ignored.');
			end
		elseif N==1
			obj.Isotope=varargin{1};
		else
			error('FluxQon:Ion:Er:WrongNargin',...
				'Incorrect number of input arguments.');
		end
	end
end

methods (Access=protected)
	function x=setIsotope(~,x)
		ME=MException('FluxQon:Ion:Er:setIsotope:InvalidInput',...
			['Input to set the isotope must be either ',...
				'0, 1 (''odd''), or 2 (''even'').']);
		if isintegerscalar(x) && any(x==[0,1,2])
			if x==0
				x=Er.randIsotope;
			end
		elseif isstringscalar(x)
			x=strcmpi(x,{'odd','even'});
			if any(x)
				x=find(x,1,'first');
			else
				throw(ME);
			end
		else
			throw(ME);
		end
	end

	function x=getElectronSpin(~,~)
		x=.5;
	end

	function x=setElectronSpin(obj,~) %#ok<STOUT>
		error('FluxQon:Ion:Er:setElectronSpin:DeniedAccess',...
			'The electron spin of %s objects is a constant.',class(obj));
	end

	function x=getNuclearSpin(obj,~)
		switch obj.Isotope
			case 1
				x=3.5;
			case 2
				x=0;
		end
	end

	function x=setNuclearSpin(obj,~) %#ok<STOUT>
		error('FluxQon:Ion:Er:setNuclearSpin:DeniedAccess',...
			'The nuclear spin of %s objects is a constant.',class(obj));
	end

	function x=getNuclearZeemanTensor(~,~)
		x=-0.1618;
	end

	function x=setNuclearZeemanTensor(obj,~) %#ok<STOUT>
		error('FluxQon:Ion:Er:setElectronSpin:DeniedAccess',...
			'The nuclear Zeeman tensor of %s objects is a constant.',class(obj));
	end
end

methods (Static=true)
	x=randIsotope(varargin)
end

end
