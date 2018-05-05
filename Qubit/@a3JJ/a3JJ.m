classdef a3JJ < FluxQubit
%% 3JJ Flux Qubit
%
% Constructor syntax:
%  obj=a3JJ() creates a a3JJ object with default properties.
%
%  obj=a3JJ(wz,wx) creates a a3JJ object with with a frequency of wz and a
%    tunneling frequency of wx.
%
%  obj=a3JJ(EC,EJ,alpha) additionally sets the alpha and beta coefficients.
%    beta is set to equal to alpha.
%
%  obj=a3JJ(EC,EJ,alpha,beta) additionally sets the alpha and beta coefficients.
%
%  obj=a3JJ('PropertyName',PropertyValue,...) sets the properties of the object
%    in Name-Value pair syntax.
%
% Requires package:
%  - MatCommon_v1.0.0+
%  - PhysConst_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: FluxQubit.
%
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 01/07/2016
% Last modified: 15/06/2017

properties
	AlphaCoefficient=.75;
	BetaCoefficient=.75;
end

properties (Dependent=true, SetAccess=private)
	CapacitancePlus
	CapacitanceMinus
	CoulombPlus
	CoulombMinus
end

methods
	function obj=a3JJ(varargin)
		if nargin==0
			% nothing is true
		elseif isstringscalar(varargin{1})
			n=1;
			while n<nargin
				mp=findprop(obj,varargin{n});
				if isempty(mp)
					error('FluxQon:Qubit:a3JJ:InvalidInput',...
						'''%s'' is not a property of the a3JJ class.',varargin{n});
				elseif strcmpi(mp.SetAccess,'public')
					obj.(varargin{n})=varargin{n+1};
				else
					error('FluxQon:Qubit:a3JJ:DeniedAccess',...
						'Denied access to set the ''%s'' property.',varargin{n});
				end
				n=n+2;
			end
			if n<=nargin
				warning('FluxQon:Qubit:a3JJ:IgnoredInput',...
					'An extra input argument was provided, but is ignored.');
			end
		elseif nargin==2
			obj.Frequency=varargin{1};
			obj.TunnelingFrequency=varargin{2};
		elseif nargin==3
			obj.CoulombEnergy=varargin{1};
			obj.JosephsonEnergy=varargin{2};
			obj.AlphaCoefficient=varargin{3};
			obj.BetaCoefficient=varargin{3};
		elseif nargin==4
			obj.CoulombEnergy=varargin{1};
			obj.JosephsonEnergy=varargin{2};
			obj.AlphaCoefficient=varargin{3};
			obj.BetaCoefficient=varargin{4};
		else
			error('FluxQon:Qubit:a3JJ:WrongNargin',...
				'Incorrect number of input arguments.');
		end
	end

	function set.AlphaCoefficient(obj,x)
		if isrealscalar(x) && x>0 && x<1
			obj.AlphaCoefficient=x;
		else
			error('FluxQon:Qubit:a3JJ:setAlphaCoefficient:InvalidInput',...
				['Input to set the alpha coefficient must be ',...
					'a real scalar between 0 and 1.']);
		end
	end

	function set.BetaCoefficient(obj,x)
		if isrealscalar(x) && x>0 && x<1
			obj.BetaCoefficient=x;
		else
			error('FluxQon:Qubit:a3JJ:setBetaCoefficient:InvalidInput',...
				['Input to set the beta coefficient must be ',...
					'a real scalar between 0 and 1.']);
		end
	end

	function x=get.CapacitancePlus(obj)
		x=(2*Constant.ElementaryCharge)^2./obj.CoulombPlus;
	end

	function x=get.CapacitanceMinus(obj)
		x=(2*Constant.ElementaryCharge)^2./obj.CoulombMinus;
	end

	function x=get.CoulombPlus(obj)
		x=obj.CoulombEnergy/2/(1+2*obj.BetaCoefficient);
	end

	function x=get.CoulombMinus(obj)
		x=obj.CoulombEnergy/2;
	end

	[E,W,dX,dY]=solveTISE(obj,N)
end

methods (Access=protected)
	function autoSetEnergy(obj)
		try
			obj.solveTISE();
		catch ME1
			if isempty(regexpi(ME1.identifier,':MissingData$','once'))
				rethrow(ME1);
			end
		end
	end

	function autoSetTunnelingEnergy(obj)
		try
			obj.solveTISE();
		catch ME1
			if isempty(regexpi(ME1.identifier,':MissingData$','once'))
				rethrow(ME1);
			end
		end
	end
end

end
