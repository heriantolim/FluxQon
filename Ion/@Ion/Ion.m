classdef Ion < handle & Object3D
%% Ion in a Microwave Field
%
% Constructor syntax:
%  obj=Ion() creates a Ion object with default properties.
%
%  obj=Ion(S) sets the electron spin to S.
%
%  obj=Ion(S,I) additionally sets the nuclear spin to I.
%
%  obj=Ion(S,gS,I,gI) additionally sets the electron (nuclear) Zeeman tensor
%    to gS (gI).
%
%  obj=Ion('PropertyName',PropertyValue,...) sets the properties of the object
%    in Name-Value pair syntax.
%
% Properties:
%  - 
%  -
%  - RoundRelTol : Relative tolerance for rounding up, specified as a positive
%                  real scalar that is less than 1. This tolerance is used in
%					    removing the trailling digits that make a Hamiltonian matrix
%					    non-Hermitian. Set to 0 to disable the rounding.
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
% First created: 15/12/2015
% Last modified: 17/06/2017

properties
	NumIons=1;
	FockDimension=3;
	MagneticField=[0;0;0];
	Temperature=0;
	RoundRelTol=1e-8;
end

properties% redefinable
	Multiplet=1;
	MultipletEnergy=0;
	Isotope=0;
	ElectronSpin=0;
	NuclearSpin=0;
	ElectronZeemanTensor=0;
	NuclearZeemanTensor=0;
	HyperfineTensor=0;
	QuadrupoleTensor=0;
	CouplingStrength=0;
	LineStrength=0;
	DecayRate=0;
end

properties (Dependent=true, SetAccess=private)
	Frequency
	NumLevels
	HilbertDimension
end

properties (SetAccess=protected)
	Energy
	Eigenstate
end

methods
	function obj=Ion(varargin)
		if nargin==0
			% nothing is true
		elseif isstringscalar(varargin{1})
			n=1;
			while n<nargin
				mp=findprop(obj,varargin{n});
				if isempty(mp)
					error('FluxQon:Ion:InvalidInput',...
						'''%s'' is not a property of the Ion class.',varargin{n});
				elseif strcmpi(mp.SetAccess,'public')
					obj.(varargin{n})=varargin{n+1};
				else
					error('FluxQon:Ion:DeniedAccess',...
						'Denied access to set the ''%s'' property.',varargin{n});
				end
				n=n+2;
			end
			if n<=nargin
				warning('FluxQon:Ion:IgnoredInput',...
					'An extra input argument was provided, but is ignored.');
			end
		elseif nargin==1
			obj.ElectronSpin=varargin{1};
		elseif nargin==2
			obj.ElectronSpin=varargin{1};
			obj.NuclearSpin=varargin{2};
		elseif nargin==4
			obj.ElectronSpin=varargin{1};
			obj.ElectronZeemanTensor=varargin{2};
			obj.NuclearSpin=varargin{3};
			obj.NuclearZeemanTensor=varargin{4};
		else
			error('FluxQon:Ion:WrongNargin',...
				'Incorrect number of input arguments.');
		end
	end

	function set.NumIons(obj,x)
		if isintegerscalar(x) && x>0
			obj.NumIons=x;
		else
			error('FluxQon:Ion:setNumIons:InvalidInput',...
				['Input to set the number of ions must be ',...
					'a positive integer scalar.']);
		end
	end

	function set.FockDimension(obj,x)
		if isintegerscalar(x) && x>0
			obj.FockDimension=x;
		else
			error('FluxQon:Ion:setNumIons:InvalidInput',...
				['Input to set the Fock dimension must be ',...
					'a positive integer scalar.']);
		end
	end

	function set.MagneticField(obj,x)
		if isrealscalar(x)
			obj.MagneticField=[0;0;x];
		elseif isrealvector(x) && numel(x)==3
			obj.MagneticField=reshape(x,3,1);
		else
			error('FluxQon:Ion:setMagneticField:InvalidInput',...
				['Input to set the magnetic field must be either ',...
					'a real scalar or a real vector of length 3.']);
		end
	end

	function set.Temperature(obj,x)
		if isrealscalar(x) && x>=0
			obj.Temperature=x;
		else
			error('FluxQon:Ion:setTemperature:InvalidInput',...
				'Input to set the temperature must be a positive real scalar.');
		end
	end

	function set.RoundRelTol(obj,x)
		if isrealscalar(x) && x>=0 && x<1
			obj.RoundRelTol=x;
		else
			error('FluxQon:Ion:setRelTol:InvalidInput',...
				['Input to set the relative error tolerance must be ',...
					'a positive real scalar less than 1.']);
		end
	end

	function x=get.Multiplet(obj)
		x=obj.Multiplet;
		x=getMultiplet(obj,x);
	end

	function set.Multiplet(obj,x)
		obj.Multiplet=setMultiplet(obj,x);
	end

	function x=get.MultipletEnergy(obj)
		x=obj.MultipletEnergy;
		x=getMultipletEnergy(obj,x);
	end

	function set.MultipletEnergy(obj,x)
		obj.MultipletEnergy=setMultipletEnergy(obj,x);
	end

	function x=get.Isotope(obj)
		x=obj.Isotope;
		x=getIsotope(obj,x);
	end

	function set.Isotope(obj,x)
		obj.Isotope=setIsotope(obj,x);
	end

	function x=get.ElectronSpin(obj)
		x=obj.ElectronSpin;
		x=getElectronSpin(obj,x);
	end

	function set.ElectronSpin(obj,x)
		obj.ElectronSpin=setElectronSpin(obj,x);
	end

	function x=get.NuclearSpin(obj)
		x=obj.NuclearSpin;
		x=getNuclearSpin(obj,x);
	end

	function set.NuclearSpin(obj,x)
		obj.NuclearSpin=setNuclearSpin(obj,x);
	end

	function x=get.ElectronZeemanTensor(obj)
		x=obj.ElectronZeemanTensor;
		x=getElectronZeemanTensor(obj,x);
	end

	function set.ElectronZeemanTensor(obj,x)
		obj.ElectronZeemanTensor=setElectronZeemanTensor(obj,x);
	end

	function x=get.NuclearZeemanTensor(obj)
		x=obj.NuclearZeemanTensor;
		x=getNuclearZeemanTensor(obj,x);
	end

	function set.NuclearZeemanTensor(obj,x)
		obj.NuclearZeemanTensor=setNuclearZeemanTensor(obj,x);
	end

	function x=get.HyperfineTensor(obj)
		x=obj.HyperfineTensor;
		x=getHyperfineTensor(obj,x);
	end

	function set.HyperfineTensor(obj,x)
		obj.HyperfineTensor=setHyperfineTensor(obj,x);
	end

	function x=get.QuadrupoleTensor(obj)
		x=obj.QuadrupoleTensor;
		x=getQuadrupoleTensor(obj,x);
	end

	function set.QuadrupoleTensor(obj,x)
		obj.QuadrupoleTensor=setQuadrupoleTensor(obj,x);
	end
	
	function x=get.CouplingStrength(obj)
		x=obj.CouplingStrength;
		x=getCouplingStrength(obj,x);
	end

	function set.CouplingStrength(obj,x)
		obj.CouplingStrength=setCouplingStrength(obj,x);
	end

	function x=get.LineStrength(obj)
		x=obj.LineStrength;
		x=getLineStrength(obj,x);
	end

	function set.LineStrength(obj,x)
		obj.LineStrength=setLineStrength(obj,x);
	end
	
	function x=get.DecayRate(obj)
		x=obj.DecayRate;
		x=getDecayRate(obj,x);
	end

	function set.DecayRate(obj,x)
		obj.DecayRate=setDecayRate(obj,x);
	end

	function x=get.Energy(obj)
		x=obj.Energy;
		if isempty(x)
			obj.Hamiltonian('S');
			x=obj.Energy;
		end
	end

	function x=get.Eigenstate(obj)
		x=obj.Eigenstate;
		if isempty(x)
			obj.Hamiltonian('S');
			x=obj.Eigenstate;
		end
	end

	function x=get.Frequency(obj)
		x=obj.Energy/Constant.ReducedPlanck;
	end

	function x=get.NumLevels(obj)
		x=(2*obj.ElectronSpin+1)*(2*obj.NuclearSpin+1)*numel(obj.Multiplet);
	end
	
	function x=get.HilbertDimension(obj)
		x=obj.NumLevels;
		if obj.NumIons>1
			x=obj.FockDimension^(x-1);
		end
	end

	A=Annihilation(obj,varargin)
	C=Creation(obj,varargin)
	N=Number(obj,varargin)
	T=Transition(obj,varargin)
	SO=ElectronSpinOperator(obj,varargin)
	IO=NuclearSpinOperator(obj,varargin)
	H=ZeemanHamiltonian(obj,varargin)
	H=HyperfineHamiltonian(obj,varargin)
	H=QuadrupoleHamiltonian(obj,varargin)
	H=Hamiltonian(obj,varargin)
	L=Lindblad(obj,varargin)
end

methods (Access=protected)
	function x=getMultiplet(~,x)
	end

	function x=setMultiplet(~,x)
		assert(isintegervector(x) && all(x>0),...
			'FluxQon:Ion:setMultiplet:InvalidInput',...
			'Input to set the multiplet must be a positive integer vector.');
		x=x(:);
	end

	function x=getMultipletEnergy(~,x)
	end

	function x=setMultipletEnergy(~,x)
		assert(isrealvector(x),...
			'FluxQon:Ion:setMultipletEnergy:InvalidInput',...
			'Input to set the multiplet energy must be a real vector.');
		x=x(:);
	end

	function x=getIsotope(~,x)
	end

	function x=setIsotope(~,x)
		assert(isintegerscalar(x) && x>=0,...
			'FluxQon:Ion:setIsotope:InvalidInput',...
			'Input to set the isotope must be a positive integer scalar.');
	end

	function x=getElectronSpin(~,x)
	end

	function x=setElectronSpin(~,x)
		assert(isnumeric(x) && isintegerscalar(2*x) && x>=0,...
			'FluxQon:setElectronSpin:InvalidInput',...
			'Input to set the electron spin must be a positive half-integer.');
	end

	function x=getNuclearSpin(~,x)
	end

	function x=setNuclearSpin(~,x)
		assert(isnumeric(x) && isintegerscalar(2*x) && x>=0,...
			'FluxQon:setNuclearSpin:InvalidInput',...
			'Input to set the nuclear spin must be a positive half-integer.');
	end

	function x=getElectronZeemanTensor(~,x)
	end

	function x=setElectronZeemanTensor(~,x)
		assert(isrealscalar(x) || (isrealmatrix(x) && all(size(x)==3)),...
			'FluxQon:setElectronZeemanTensor:InvalidInput',...
			['Input to set the electron Zeeman tensor must be either ',...
				'a real scalar or a real matrix of size 3x3.']);
		if all(x(:)==0)
			x=0;
		end
	end

	function x=getNuclearZeemanTensor(~,x)
	end

	function x=setNuclearZeemanTensor(~,x)
		assert(isrealscalar(x) || (isrealmatrix(x) && all(size(x)==3)),...
			'FluxQon:setNuclearZeemanTensor:InvalidInput',...
			['Input to set the nuclear Zeeman tensor must be either ',...
				'a real scalar or a real matrix of size 3x3.']);
		if all(x(:)==0)
			x=0;
		end
	end

	function x=getHyperfineTensor(~,x)
	end

	function x=setHyperfineTensor(~,x)
		assert(isrealscalar(x) || (isrealmatrix(x) && all(size(x)==3)),...
			'FluxQon:setHyperfineTensor:InvalidInput',...
			['Input to set the hyperfine tensor must be either ',...
				'a real scalar or a real matrix of size 3x3.']);
		if all(x(:)==0)
			x=0;
		end
	end

	function x=getQuadrupoleTensor(~,x)
	end

	function x=setQuadrupoleTensor(~,x)
		assert(isrealscalar(x) || (isrealmatrix(x) && all(size(x)==3)),...
			'FluxQon:setQuadrupoleTensor:InvalidInput',...
			['Input to set the quadrupole tensor must be either ',...
				'a real scalar or a real matrix of size 3x3.']);
		if all(x(:)==0)
			x=0;
		end
	end
	
	function x=getCouplingStrength(~,x)
	end
	
	function x=setCouplingStrength(~,x)
		if iscomplexvector(x)
			x=vec2trlherm(x);
		elseif isempty(x) || ~ishermitian(x)
			error('FluxQon:Ion:setCouplingStrength:InvalidInput',...
				['Input to set the coupling strength must be either ',...
					'a complex vector or a hermitian matrix.']);
		end
	end

	function x=getLineStrength(~,x)
	end

	function x=setLineStrength(~,x)
		if iscomplexvector(x)
			x=vec2trlherm(x);
		elseif isempty(x) || ~ishermitian(x)
			error('FluxQon:Ion:setLineStrength:InvalidInput',...
				['Input to set the line strength must be either ',...
					'a complex vector or a hermitian matrix.']);
		end
	end
	
	function x=getDecayRate(~,x)
	end
	
	function x=setDecayRate(~,x)
		assert(~isempty(x) && isrealmatrix(x) && all(x(:)>=0) ...
				&& (isvector(x) || issymmetric(x)),...
			'FluxQon:Ion:setDecayRate:InvalidInput',...
			['Input to set the decay rate must be either ',...
				'a positive real vector or a positive real symmetric matrix.']);
	end

	T=rotateTensor(obj,T)
end

end
