classdef FluxQubit < Qubit & Object2D
%% Flux Qubit
%
% Constructor syntax:
%  obj=FluxQubit() creates a FluxQubit object with default properties.
%
%  obj=FluxQubit(wz,wx) creates a FluxQubit object with a frequency of wz and a
%    tunneling frequency of wx.
%
%  obj=FluxQubit(EC,EL,EJ) creates a FluxQubit object with a Coulomb energy EC,
%    an inductance energy EL, and a Josephson energy EJ.
%
%  obj=FluxQubit('PropertyName',PropertyValue,...) sets the properties of the
%    object in Name-Value pair syntax.
%
% Requires package:
%  - Common_v1.0.0+
%  - PhysConst_v1.0.0+
%  - QuantMech_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: Qubit.
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 15/06/2017
% Last modified: 15/06/2017

properties
	CoulombEnergy
	InductanceEnergy
	JosephsonEnergy
	PhaseBias=pi;
	PhaseTune=0;
end

properties (Dependent=true)
	Capacitance
	SelfInductance
	CriticalCurrent
	FluxBias
	FluxTune
end

methods
	function obj=FluxQubit(varargin)
		if nargin==0
			% nothing is true
		elseif isstringscalar(varargin{1})
			n=1;
			while n<nargin
				mp=findprop(obj,varargin{n});
				if isempty(mp)
					error('FluxQon:Qubit:FluxQubit:InvalidInput',...
						'''%s'' is not a property of the FluxQubit class.',...
						varargin{n});
				elseif strcmpi(mp.SetAccess,'public')
					obj.(varargin{n})=varargin{n+1};
				else
					error('FluxQon:Qubit:FluxQubit:DeniedAccess',...
						'Denied access to set the ''%s'' property.',varargin{n});
				end
				n=n+2;
			end
			if n<=nargin
				warning('FluxQon:Qubit:FluxQubit:IgnoredInput',...
					'An extra input argument was provided, but is ignored.');
			end
		elseif nargin==2
			obj.Frequency=varargin{1};
			obj.TunnelingFrequency=varargin{2};
		elseif nargin==3
			obj.CoulombEnergy=varargin{1};
			obj.InductanceEnergy=varargin{2};
			obj.JosephsonEnergy=varargin{3};
		else
			error('FluxQon:Qubit:FluxQubit:WrongNargin',...
				'Incorrect number of input arguments.');
		end
	end

	function set.CoulombEnergy(obj,x)
		if isrealscalar(x) && x>=0
			obj.CoulombEnergy=x;
		else
			error('FluxQon:Qubit:FluxQubit:setCoulombEnergy:InvalidInput',...
				'Input to set the Coulomb energy must be a positive real scalar.');
		end
	end

	function set.InductanceEnergy(obj,x)
		if isrealscalar(x) && x>=0
			obj.InductanceEnergy=x;
		else
			error('FluxQon:Qubit:FluxQubit:setInductanceEnergy:InvalidInput',...
				['Input to set the Inductance energy must be ',...
					'a positive real scalar.']);
		end
	end

	function set.JosephsonEnergy(obj,x)
		if isrealscalar(x) && x>=0
			obj.JosephsonEnergy=x;
		else
			error('FluxQon:Qubit:FluxQubit:setJosephsonEnergy:InvalidInput',...
				['Input to set the Josephson energy must be ',...
					'a positive real scalar.']);
		end
	end

	function set.PhaseBias(obj,x)
		if isrealscalar(x)
			obj.PhaseBias=x;
		else
			error('FluxQon:Qubit:FluxQubit:setPhaseBias:InvalidInput',...
				'Input to set the phase bias must be a real scalar.');
		end
	end

	function set.PhaseTune(obj,x)
		if isrealscalar(x)
			obj.PhaseTune=x;
		else
			error('FluxQon:Qubit:FluxQubit:setPhaseTune:InvalidInput',...
				'Input to set the phase tune must be a real scalar.');
		end
	end

	function x=get.Capacitance(obj)
		x=(2*Constant.ElementaryCharge)^2./obj.CoulombEnergy;
	end

	function set.Capacitance(obj,x)
		if isrealscalar(x) && x>=0
			obj.CoulombEnergy=(2*Constant.ElementaryCharge)^2/x;
		else
			error('FluxQon:Qubit:FluxQubit:setCapacitance:InvalidInput',...
				'Input to set the capacitance must be a positive real scalar.');
		end
	end

	function x=get.SelfInductance(obj)
		x=(Constant.FluxQuantum/2/pi)^2./obj.InductanceEnergy;
	end

	function set.SelfInductance(obj,x)
		if isrealscalar(x) && x>=0
			obj.InductanceEnergy=(Constant.FluxQuantum/2/pi)^2/x;
		else
			error('FluxQon:Qubit:FluxQubit:setSelfInductance:InvalidInput',...
				'Input to set the self inductance must be a positive real scalar.');
		end
	end

	function x=get.CriticalCurrent(obj)
		x=2*pi/Constant.FluxQuantum*obj.JosephsonEnergy;
	end

	function set.CriticalCurrent(obj,x)
		if isrealscalar(x) && x>=0
			obj.JosephsonEnergy=Constant.FluxQuantum/2/pi*x;
		else
			error('FluxQon:Qubit:FluxQubit:setCriticalCurrent:InvalidInput',...
				['Input to set the critical current must be a positive real ',...
					'scalar.']);
		end
	end

	function x=get.FluxBias(obj)
		x=Constant.FluxQuantum/2/pi.*obj.PhaseBias;
	end

	function set.FluxBias(obj,x)
		if isrealscalar(x)
			obj.PhaseBias=2*pi*x/Constant.FluxQuantum;
		else
			error('FluxQon:Qubit:FluxQubit:setFluxBias:InvalidInput',...
				'Input to set the flux bias must be real scalar.');
		end
	end

	function x=get.FluxTune(obj)
		x=Constant.FluxQuantum/2/pi.*obj.PhaseTune;
	end

	function set.FluxTune(obj,x)
		if isrealscalar(x)
			obj.PhaseTune=2*pi*x/Constant.FluxQuantum;
		else
			error('FluxQon:Qubit:FluxQubit:setFluxTune:InvalidInput',...
				'Input to set the flux tune must be real scalar.');
		end
	end
end

end
