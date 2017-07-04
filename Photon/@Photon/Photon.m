classdef Photon < handle
%% Photon
%
% Constructor syntax:
%  obj=Photon() creates a Photon object with a random angular frequency.
%
%  obj=Photon(w) creates a Photon object with an angular frequency w.
%
%  obj=Photon(w,F) additionally sets the dimension of the Fock space to F.
%
%  obj=Photon(w,F,V) additionally sets the confinement volume to V.
%
%  obj=Photon('PropertyName',PropertyValue,...) sets the properties of the
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
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 10/06/2017
% Last modified: 15/06/2017

properties
	FockDimension=3;
	Energy
	PumpEnergy=0;
	PumpPower=0;
	Temperature=0;
	DecayRate=0;
	IncoherentPumpRate=0;
	ConfinementVolume=1e-12;
end

properties (Dependent=true)
	Frequency% Angular frequency
	Wavelength
	PumpFrequency% Angular frequency
	PumpRate
	ElectricAmplitude
	MagneticAmplitude
end

properties (Dependent=true, SetAccess=private)
	HilbertDimension
end

methods
	function obj=Photon(varargin)
		if nargin==0
			R=obj.FrequencyRange;
			obj.Frequency=R(1)+(R(2)-R(1))^rand;
		elseif isstringscalar(varargin{1})
			n=1;
			while n<nargin
			mp=findprop(obj,varargin{n});
				if isempty(mp)
					error('FluxQon:Photon:InvalidInput',...
						'''%s'' is not a property of the Photon class.',varargin{n});
				elseif strcmpi(mp.SetAccess,'public')
					obj.(varargin{n})=varargin{n+1};
				else
					error('FluxQon:Photon:DeniedAccess',...
						'Denied access to set the ''%s'' property.',varargin{n});
				end
				n=n+2;
			end
			if n<=nargin
				warning('FluxQon:Photon:IgnoredInput',...
					'An extra input argument was provided, but is ignored.');
			end
		elseif nargin==1
			obj.Frequency=varargin{1};
		elseif nargin==2
			obj.Frequency=varargin{1};
			obj.FockDimension=varargin{2};
		elseif nargin==3
			obj.Frequency=varargin{1};
			obj.FockDimension=varargin{2};
			obj.ConfinementVolume=varargin{3};
		else
			error('FluxQon:Photon:WrongNargin',...
				'Incorrect number of input arguments.');
		end
	end

	function set.FockDimension(obj,x)
		if isintegerscalar(x) && x>1
			obj.FockDimension=x;
		else
			error('FluxQon:Photon:setFockDimension:InvalidInput',...
				['Input to set the Fock dimension must be an integer scalar ',...
					'greater than 1.']);
		end
	end

	function set.Energy(obj,x)
		R=Constant.ReducedPlanck*obj.FrequencyRange;
		if isrealscalar(x) && x>R(1) && x<R(2)
			obj.Energy=x;
		else
			error('FluxQon:Photon:setEnergy:InvalidInput',...
				['Input to set the energy must be a real scalar ',...
					'between %g and %g.'],R(1),R(2));
		end
	end
	
	function set.PumpEnergy(obj,x)
		if isrealscalar(x) && x>=0
			obj.PumpEnergy=x;
		else
			error('FluxQon:Photon:setPumpEnergy:InvalidInput',...
				'Input to set the pump energy must be a positive real scalar.');
		end
	end
	
	function set.PumpPower(obj,x)
		if isrealscalar(x) && x>=0
			obj.PumpPower=x;
		else
			error('FluxQon:Photon:setPumpPower:InvalidInput',...
				'Input to set the pump power must be a positive real scalar.');
		end
	end
	
	function set.Temperature(obj,x)
		if isrealscalar(x) && x>=0
			obj.Temperature=x;
		else
			error('FluxQon:Photon:setTemperature:InvalidInput',...
				'Input to set the temperature must be a positive real scalar.');
		end
	end
	
	function set.DecayRate(obj,x)
		if isrealscalar(x) && x>=0
			obj.DecayRate=x;
		else
			error('FluxQon:Photon:setDecayRate:InvalidInput',...
				'Input to set the decay rate must be a positive real scalar.');
		end
	end
	
	function set.IncoherentPumpRate(obj,x)
		if isrealscalar(x) && x>=0
			obj.IncoherentPumpRate=x;
		else
			error('FluxQon:Photon:setIncoherentPumpRate:InvalidInput',...
				['Input to set the incoherent pump rate must be ',...
					'a positive real scalar.']);
		end
	end

	function set.ConfinementVolume(obj,x)
		if isrealscalar(x) && x>0
			obj.ConfinementVolume=x;
		else
			error('FluxQon:Photon:setConfinementVolume:InvalidInput',...
				['Input to set the confinement volume must be ',...
					'a positive real scalar.']);
		end
	end

	function x=get.Frequency(obj)
		x=obj.Energy/Constant.ReducedPlanck;
	end

	function set.Frequency(obj,x)
		R=obj.FrequencyRange;
		if isrealscalar(x) && x>R(1) && x<R(2)
			obj.Energy=Constant.ReducedPlanck*x;
		else
			error('FluxQon:Photon:setFrequency:InvalidInput',...
				['Input to set the frequency must be a real scalar ',...
					'between %g and %g.'],R(1),R(2));
		end
	end

	function x=get.Wavelength(obj)
		x=2*pi*Constant.LightSpeed/obj.Frequency;
	end

	function set.Wavelength(obj,x)
		R=2*pi*Constant.LightSpeed./obj.FrequencyRange;
		if isrealscalar(x) && x>R(2) && x<R(1)
			obj.Energy=Constant.Planck*Constant.LightSpeed/x;
		else
			error('FluxQon:Photon:setWavelength:InvalidInput',...
				['Input to set the wavelength must be a real scalar ',...
					'between %g and %g.'],R(2),R(1));
		end
	end
	
	function x=get.PumpFrequency(obj)
		x=obj.PumpEnergy/Constant.ReducedPlanck;
	end
	
	function set.PumpFrequency(obj,x)
		if isrealscalar(x) && x>=0
			obj.PumpEnergy=Constant.ReducedPlanck*x;
		else
			error('FluxQon:setPumpFrequency:InvalidInput',...
				'Input to set the pump frequency must be a positive real scalar.');
		end
	end
	
	function x=get.PumpRate(obj)
		x=obj.PumpPower/Constant.ReducedPlanck;
	end
	
	function set.PumpRate(obj,x)
		if isrealscalar(x) && x>=0
			obj.PumpPower=Constant.ReducedPlanck*x;
		else
			error('FluxQon:setPumpRate:InvalidInput',...
				'Input to set the pump rate must be a positive real scalar.');
		end
	end

	function x=get.ElectricAmplitude(obj)
		x=sqrt(Constant.ReducedPlanck/2/Constant.VacuumPermittivity ...
			*obj.Frequency/obj.ConfinementVolume);
	end

	function set.ElectricAmplitude(obj,x)
		if isrealscalar(x) && x>0
			obj.ConfinementVolume=Constant.ReducedPlanck/2 ...
				/Constant.VacuumPermittivity*obj.Frequency/x^2;
		else
			error('FluxQon:Photon:setElectricAmplitude:InvalidInput',...
				['Input to set the electric amplitude must be ',...
					'a positive real scalar.']);
		end
	end

	function x=get.MagneticAmplitude(obj)
		x=sqrt(Constant.ReducedPlanck/2/Constant.VacuumPermittivity ...
			/obj.Frequency/obj.ConfinementVolume);
	end

	function set.MagneticAmplitude(obj,x)
		if isrealscalar(x) && x>0
			obj.ConfinementVolume=Constant.ReducedPlanck/2 ...
				/Constant.VacuumPermittivity/obj.Frequency/x^2;
		else
			error('FluxQon:Photon:setMagneticAmplitude:InvalidInput',...
				['Input to set the magnetic amplitude must be ',...
					'a positive real scalar.']);
		end
	end

	function x=get.HilbertDimension(obj)
		x=obj.FockDimension;
	end

	A=Annihilation(obj)
	C=Creation(obj)
	N=Number(obj)
	H=Hamiltonian(obj)
	L=Lindblad(obj,varargin)
end

methods (Static=true)
	function x=FrequencyRange()
		x=[0,1.1655e44];
	end
end

end
