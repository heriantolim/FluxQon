classdef Qubit < handle
%% Qubit
%
% Constructor syntax:
%  obj=Qubit() creates a Qubit object with default properties.
%
%  obj=Qubit(wz) creates a Qubit object with a frequency of wz.
%
%  obj=Qubit(wz,wx) creates a Qubit object with a frequency of wz and a
%    tunneling frequency of wx.
%
%  obj=Qubit('PropertyName',PropertyValue,...) sets the properties of the object
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
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 18/12/2015
% Last modified: 17/06/2017

properties
	Energy
	TunnelingEnergy
	Temperature=0;
	DecayRate=0;
	DephasingRate=0;
end

properties (Dependent=true)
	Frequency
	TunnelingFrequency
end

properties (Constant=true)
	HilbertDimension=2;
end

methods
	function obj=Qubit(varargin)
		if nargin==0
			% nothing is true
		elseif isstringscalar(varargin{1})
			n=1;
			while n<nargin
				mp=findprop(obj,varargin{n});
				if isempty(mp)
					error('FluxQon:Qubit:InvalidInput',...
						'''%s'' is not a property of the Qubit class.',varargin{n});
				elseif strcmpi(mp.SetAccess,'public')
					obj.(varargin{n})=varargin{n+1};
				else
					error('FluxQon:Qubit:DeniedAccess',...
						'Denied access to set the ''%s'' property.',varargin{n});
				end
				n=n+2;
			end
			if n<=nargin
				warning('FluxQon:Qubit:IgnoredInput',...
					'An extra input argument was provided, but is ignored.');
			end
		elseif nargin==1
			obj.Frequency=varargin{1};
		elseif nargin==2
			obj.Frequency=varargin{1};
			obj.TunnelingFrequency=varargin{2};
		else
			error('FluxQon:Qubit:WrongNargin',...
				'Incorrect number of input arguments.');
		end
	end

	function x=get.Energy(obj)
		x=obj.Energy;
		if isempty(x)
			obj.autoSetEnergy;
			x=obj.Energy;
		end
	end

	function set.Energy(obj,x)
		if isrealscalar(x) && x>=0
			obj.Energy=x;
		else
			error('FluxQon:Qubit:setEnergy:InvalidInput',...
				'Input to set the energy must be a positive real scalar.');
		end
	end

	function x=get.TunnelingEnergy(obj)
		x=obj.TunnelingEnergy;
		if isempty(x)
			obj.autoSetTunnelingEnergy;
			x=obj.TunnelingEnergy;
		end
	end

	function set.TunnelingEnergy(obj,x)
		if isrealscalar(x)
			obj.TunnelingEnergy=x;
		else
			error('FluxQon:Qubit:setTunnelingEnergy:InvalidInput',...
				'Input to set the tunneling energy must be a real scalar.');
		end
	end

	function set.Temperature(obj,x)
		if isrealscalar(x) && x>=0
			obj.Temperature=x;
		else
			error('FluxQon:Qubit:setTemperature:InvalidInput',...
				'Input to set the temperature must be a positive real scalar.');
		end
	end

	function set.DecayRate(obj,x)
		if isrealscalar(x) && x>=0
			obj.DecayRate=x;
		else
			error('FluxQon:Qubit:setDecayRate:InvalidInput',...
				'Input to set the decay rate must be a positive real scalar.');
		end
	end

	function set.DephasingRate(obj,x)
		if isrealscalar(x) && x>=0
			obj.DephasingRate=x;
		else
			error('FluxQon:Qubit:setDephasingRate:InvalidInput',...
				'Input to set the dephasing rate must be a positive real scalar.');
		end
	end

	function x=get.Frequency(obj)
		x=obj.Energy/Constant.ReducedPlanck;
	end

	function set.Frequency(obj,x)
		if isrealscalar(x) && x>=0
			obj.Energy=Constant.ReducedPlanck*x;
		else
			error('FluxQon:Qubit:setFrequency:InvalidInput',...
				'Input to set the frequency must be a positive real scalar.');
		end
	end

	function x=get.TunnelingFrequency(obj)
		x=obj.TunnelingEnergy/Constant.ReducedPlanck;
	end

	function set.TunnelingFrequency(obj,x)
		if isrealscalar(x)
			obj.TunnelingEnergy=Constant.ReducedPlanck*x;
		else
			error('FluxQon:Qubit:setTunnelingFrequency:InvalidInput',...
				'Input to set the tunneling frequency must be a real scalar.');
		end
	end

	H=Hamiltonian(obj)
	L=Lindblad(obj,varargin)
end

methods (Static=true)
	A=Annihilation()
	C=Creation()
	N=Number()
end

methods (Access=protected)
	function autoSetEnergy(obj) %#ok<MANU>
	end

	function autoSetTunnelingEnergy(obj) %#ok<MANU>
	end
end

end
