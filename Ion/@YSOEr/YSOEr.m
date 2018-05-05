classdef YSOEr < YSO & Er
%% Y2SiO5:Er3+ in a Microwave Field
%
% Constructor syntax:
%  obj=YSOEr(), obj=YSOEr(0) creates an YSOEr object with a random isotope,
%  site, and class. The isotope is randomed according to its natural abundance:
%  22.95% chance being odd, and the rest being even. The site and class have a
%  50:50 chance of being 1 or 2.
%
%  obj=Er(1), obj=Er('odd') creates an YSOEr object with an odd isotope, a
%  random site and a random class.
%
%  obj=Er(2), obj=Er('even') creates an YSOEr object with an even isotope, a
%  random site and a random class.
%
%  obj=YSOEr([a,b,c])
%
%  obj=YSOEr(isotope,site)
%  obj=YSOEr(isotope,[a,b,c])
%
%  obj=YSOEr(isotope,site,class)
%  obj=YSOEr(isotope,site,[a,b,c])
%
%  obj=YSOEr(isotope,site,class,B), where B is scalar.
%  obj=YSOEr(isotope,site,class,[a,b,c])
%
%  obj=YSOEr(isotope,site,class,[a,b,c],B), where B can be scalar or vector.
%
%  obj=YSOEr(N,F,...), where N>2, creates an YSOEr object representing a
%  collection of N YSOEr ions. F is the dimension of the annihilation-and-
%  creation operators used to approximate the spin operators.
%
%  obj=YSOEr('PropertyName',PropertyValue,...) sets the properties of the object
%  in Name-Value pair syntax.
%
% Input arguments:
%  isotope : Er3+ isotope, specified as an integer of either 1 or 2. Can also
%            be set as either 'odd' (1) or 'even' (2). Defaults to random with a
%            22.95% chance being 1. Set 0 to random.
%
%  site    : YSO site, specified as an integer of either 1 or 2. Defaults to
%            random with a 50% chance being 1. Set 0 to random.
%
%  class   : Class of the magnetic orientation, specified as an integer of
%            either 1 or 2. Defaults to random with a 50% chance being 1. Set 0
%            to random.
%
%  [a,b,c] : Rotation of the YSO in the 3D space. a, b, and c are angles in
%            radians.
%
%  B       : External magnetic field, specified as a real scalar or a real
%            vector of length 3. When specified as a scalar, B is set to be
%            equal to [0;0;B].
%
% Requires package:
%  - MatCommon_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: Ion, Er.
%
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 05/12/2015
% Last modified: 17/06/2017

properties
	Site=0;% either 1 or 2. Set 0 to random. (50%,50%)
	Class=0;% either 1 or 2. Set 0 to random. (50%,50%)
end

methods
	function obj=YSOEr(varargin)
		N=nargin;
		if N>1 && isintegerscalar(varargin{1}) && varargin{1}>2
			obj.NumIons=varargin{1};
			obj.FockDimension=varargin{2};
			varargin=varargin(3:N);
			N=N-2;
		end
		if N==0
			% nothing is true
		elseif isstringscalar(varargin{1})
			n=1;
			while n<N
				mp=findprop(obj,varargin{n});
				if isempty(mp)
					error('FluxQon:Ion:YSOEr:InvalidInput',...
						'''%s'' is not a property of the YSOEr class.',varargin{n});
				elseif strcmpi(mp.SetAccess,'public')
					obj.(varargin{n})=varargin{n+1};
				else
					error('FluxQon:Ion:YSOEr:DeniedAccess',...
						'Denied access to set the ''%s'' property.',varargin{n});
				end
				n=n+2;
			end
			if n<=N
				warning('FluxQon:Ion:YSOEr:IgnoredInput',...
					'An extra input argument was provided, but is ignored.');
			end
		elseif N==1
			if isscalar(varargin{1})
				obj.Isotope=varargin{1};
			else
				obj.Rotation=varargin{1};
			end
		elseif N==2
			obj.Isotope=varargin{1};
			if isscalar(varargin{2})
				obj.Site=varargin{2};
			else
				obj.Rotation=varargin{2};
			end
		elseif N==3
			obj.Isotope=varargin{1};
			obj.Site=varargin{2};
			if isscalar(varargin{3})
				obj.Class=varargin{3};
			else
				obj.Rotation=varargin{3};
			end
		elseif N==4
			obj.Isotope=varargin{1};
			obj.Site=varargin{2};
			obj.Class=varargin{3};
			if isscalar(varargin{4})
				obj.MagneticField=varargin{4};
			else
				obj.Rotation=varargin{4};
			end
		elseif N==5
			obj.Isotope=varargin{1};
			obj.Site=varargin{2};
			obj.Class=varargin{3};
			obj.Rotation=varargin{4};
			obj.MagneticField=varargin{5};
		else
			error('FluxQon:Ion:YSOEr:WrongNargin',...
				'Incorrect number of input arguments.');
		end
		if obj.Isotope==0
			obj.Isotope=0;
		end
		if obj.Site==0
			obj.Site=0;
		end
		if obj.Class==0
			obj.Class=0;
		end
	end

	function set.Site(obj,x)
		if isintegerscalar(x) && any(x==[0,1,2])
			if x==0
				obj.Site=YSOEr.randSite;
			else
				obj.Site=x;
			end
		else
			error('FluxQon:Ion:YSOEr:setSite:InvalidInput',...
				'Input to set the site must be either 0, 1, or 2.');
		end
	end

	function set.Class(obj,x)
		if isintegerscalar(x) && any(x==[0,1,2])
			if x==0
				obj.Class=YSOEr.randClass;
			else
				obj.Class=x;
			end
		else
			error('FluxQon:Ion:YSOEr:setClass:InvalidInput',...
				'Input to set the class must be either 0, 1, or 2.');
		end
	end
end

methods (Access=protected)
	function x=setMultiplet(~,x)
		if isintegervector(x) && numel(x)<3 && all(ismember(x,[1,2]))
			x=unique(x);
		else
			error('FluxQon:Ion:YSOEr:setMultiplet:InvalidInput',...
				'Input to set the multiplet must be either 1, 2, or [1,2]');
		end
	end

	function x=getMultipletEnergy(obj,~)
		MV=obj.Multiplet;
		M=numel(MV);
		x=zeros(M,1);
		for Mi=1:M
			if MV(Mi)==2
				switch obj.Site
					case 1
						x(Mi)=1.292857e-19;
					case 2
						x(Mi)=1.290819e-19;
				end
			end
		end
	end

	function x=setMultipletEnergy(obj,~) %#ok<STOUT>
		error('FluxQon:Ion:YSOEr:setMultiPletEnergy:DeniedAccess',...
			'The multiplet energy of %s objects is a constant.',class(obj));
	end

	function x=getElectronZeemanTensor(obj,~)
		switch obj.Site
			case 1
				x={[  3.070, -3.124,  3.396;
					  -3.124,  8.156, -5.756;
					   3.396, -5.756,  5.787],...
				   [  1.950, -2.212,  3.584;
					  -2.212,  4.232, -4.986;
					   3.584, -4.986,  7.888]};
			case 2
				x={[ 14.651, -2.115,  2.552;
					  -2.115,  1.965, -0.550;
					   2.552, -0.550,  0.902],...
				   [ 12.032, -0.582,  4.518;
					  -0.582,  0.212, -0.296;
					   4.518, -0.296,  1.771]};
		end
		M=obj.Multiplet;
		x=x(M);
		M=numel(M);
		if obj.Class==2
			for Mi=1:M
				x{Mi}=diag([-1,-1,1])*x{Mi}*diag([-1,-1,1]);
			end
		end
		for Mi=1:M
			x{Mi}=obj.rotateTensor(x{Mi});
		end
		x=cat(3,x{:});
	end

	function x=setElectronZeemanTensor(obj,~) %#ok<STOUT>
		error('FluxQon:Ion:YSOEr:setElectronZeemanTensor:DeniedAccess',...
			'The electron Zeeman tensor of %s objects is a constant.',class(obj));
	end

	function x=getHyperfineTensor(obj,~)
		M=obj.Multiplet;
		if obj.Isotope==1 && ~isequal(M,2)
			switch obj.Site
				case 1
					x=[   69.35, -580.73,  248.83;
						 -580.73,  696.30, -682.49;
						  248.83, -682.49,  495.54];
				case 2
					x=[-1521.4 ,  178.11,  141.76;
						  178.11,  172.09, -212.54;
						  141.76, -212.54,  199.01];
			end
			if obj.Class==2
				x=diag([-1,-1,1])*x*diag([-1,-1,1]);
			end
			x=1e6*Constant.Planck*obj.rotateTensor(x);
			if numel(M)==2
				x=cat(3,x,zeros(3));
			end
		else
			x=0;
		end
	end

	function x=setHyperfineTensor(obj,~) %#ok<STOUT>
		error('FluxQon:Er:setHyperfineTensor:DeniedAccess',...
			'The hyperfine tensor of %s objects is a constant.',class(obj));
	end

	function x=getQuadrupoleTensor(obj,~)
		M=obj.Multiplet;
		if obj.Isotope==1 && ~isequal(M,2)
			switch obj.Site
				case 1
					x=[  21.4, -8.18, 15.27;
						 -8.18,  3.8 , -0.6 ;
						 15.27, -0.6 ,-25.2 ];
				case 2
					x=[ -3.5 ,-19.84,-24.22;
						-19.84, 50.4 , -6.73;
						-24.22, -6.73,-46.9 ];
			end
			if obj.Class==2
				x=diag([-1,-1,1])*x*diag([-1,-1,1]);
			end
			x=1e6*Constant.Planck*obj.rotateTensor(x);
			if numel(M)==2
				x=cat(3,x,zeros(3));
			end
		else
			x=0;
		end
	end

	function x=setQuadrupoleTensor(obj,~) %#ok<STOUT>
		error('FluxQon:Er:setQuadrupoleTensor:DeniedAccess',...
			'The quadrupole tensor of %s objects is a constant.',class(obj));
	end
end

methods (Static=true)
	x=randSite(varargin)
	x=randClass(varargin)
end

end
