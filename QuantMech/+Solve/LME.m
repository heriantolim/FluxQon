function [t,rho,varargout]=LME(tspan,rho0,L,varargin)
%% Solve the Lindblad Master Equation
%  [t,rho]=Solve.LME(tspan,rho0,L) solves the evolution of a quantum state as
%    governed by the LME with a constant Lindbladian superoperator L. The state
%    is evolved within the time interval tspan from the initial density matrix
%    rho0. This function utilizes the MATLAB function ode45, which is a
%    single-step solver based on the adaptive Runge-Kutta method. Use the
%    Lindblad function to form the superoperator L.
%
%	[t,rho]=Solve.LME(tspan,rho0,L,H1) takes a time-dependent Hamiltonian H1 into
%	  the equation. H1 is to be specified as a function handle that accepts the
%	  time as its first argument and optionally the density matrix as its second
%	  argument.
%
%  [t,rho]=Solve.LME(..., Name-Value Args) accepts inputs to set the options for
%    ODE solver. The input arguments are passed in the Name-Value format. See
%    the available options below.
%
%  [t,rho,te,ye,ie]=Solve.LME(...,'Events',fun) checks the solution at each
%    iteration step with the 'event' function. The event function can be used to
%    mark the specific times at which certain conditions are met and to stop the
%    iteration if necessary. Refer to the MATLAB Help page titled 'ODE Event
%    Location' for more information.
%
% Outputs:
%  t: Evaluation times, returned as a row vector. If tspan contains two
%     elements, [t0,tf], then t contains the internal evaluation points used to
%     perform the integration. Otherwise t is the same as tspan.
%
%  rho: Solutions, returned as a complex three-dimensional array. Each matrix in
%       rho( , ,i) corresponds to the solution at time t(i).
%
%  te: Time of events, returned as a row vector. The event times in te
%      correspond to the solutions returned in ye, and ie specifies which event
%      occured.
%
%  ye: Solution at time of events, returned as a complex matrix. Each column in
%      ye corresponds to the solution at the value returned in the corresponding
%      column of te.
%
%  ie: Index of vanishing event function, returned as a column vector. This
%      index tells which event occured.
%
% Required Inputs:
%  tspan: Time span of the evolution, specified as a real vector. At minimum,
%         tspan must be a two element vector [t0,tf] specifying the initial and
%         final times. To obtain solutions at specific times between t0 and tf,
%         use a longer vector of the form [t0,t1,t2,...,tf].
%
%  rho0 : Initial density matrix of the system, specified as a complex matix of
%         size DxD. The trace of rho0 must equal to 1.
%
%  L    : Lindblad superoperator, specified as a complex matrix of size D^2xD^2.
%
% Optional Inputs:
%  H1: Time-dependent Hamiltonian, specified as a function handle. The function
%      can take one or two arguments, the first is for the time, and the second
%      is for the density matrix at the given time. The function must return a
%      Hermitian matrix of size DxD.
%
% Name-Value Inputs:
%  RelTol: Relative error tolerance, specified as a positive real scalar.
%          Defaults to 1e-3. This tolerance sets the upper bound for the error
%          relative to the magnitude of each solution component. At each
%          iteration steps, the solver checks if the relative error is below
%          RelTol. If it is above, then the computation is repeated with smaller
%          time step.
%
%  AbsTol: Absolute error tolerance, specified as a positive real scalar.
%          Defaults to 1e-8. The solver also makes sure that the absolute error
%          of the solution is below this value.
%
%  InitialStep: Suggested initial step size as a fraction of tf-t0, specified as
%               a real scalar between 0 and 1. This option sets an upper
%               bound for the magnitude of the first step size that the solver
%               tries. If not specified, then the intial step is based on the
%               slope of the solution at the initial time point.
%
%  MaxStep: Maximum step size as a fraction of tf-t0, specified as a real scalar
%           between 0 and 1. Defaults to .01. This option sets an upper bound on
%           the size of any step taken by the solver. If the ODE has periodic
%           behavior, for example, then setting MaxStep to a fraction of the
%           period ensures that the solver does not enlarge the step so much
%           that it steps over an area of interest.
%
%  Events: Event function, specified as a function handle. The event function
%          specified by the handle must have the general form as outlined in the
%          MATLAB Help page titled 'ODE Event Location'.
%
%  PlotObs: Observable(s) to be plotted, specified as a complex matrix of size
%           DxD or as a cell containing DxD matrices. The latter input format
%           allows several observables to be plotted simultaneously. The solver
%           plots the expectation value of the observables at each successful
%           time step. Specifying this option slows down the computation.
%
% Requires package:
%  - Common_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: Lindblad, TDSE.
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 17/05/2017
% Last modified: 15/06/2017

%% Defaults
H1=[];
RelTol=1e-3;
AbsTol=1e-6;
InitialStep=0;
MaxStep=.01;
Events=[];
PlotObs=[];

%% Input Parsing and Validation
assert(isrealvector(tspan),...
	'QuantMech:Solve:LME:InvalidInput',...
	'Input to the time span must be a real vector.');
tspan=sort(tspan);

D=size(rho0,1);
assert(iscomplexmatrix(rho0) && size(rho0,2)==D,...
	'QuantMech:Solve:LME:InvalidInput',...
	'Input to the initial density matrix must be a complex square matrix.');

assert(iscomplexmatrix(L) && all(size(L)==D^2),...
	'QuantMech:Solve:LME:InvalidInput',...
	['Input to the Lindblad superoperator must be a complex matrix and '...
		'the matrix dimension must equal to the squared dimension of the ',...
		'density matrix.']);

N=nargin-3;
if N>0
	n=1;
	x=varargin{n};
	if isa(x,'function_handle')
		H1f=nargin(x);
		switch H1f
			case 1
				H1=x(tspan(1));
			case 2
				H1=x(tspan(1),normc(rand(D,1)));
			otherwise
				error('QuantMech:Solve:LME:InvalidInput',...
					['The function for the time-dependent Hamiltonian must take ',...
						'one or two input arguments']);
		end
		H1f=H1f==1;
		assert(ishermitian(H1) && size(H1,1)==D,...
			'QuantMech:Solve:LME:InvalidInput',...
			['The function for the time-dependent Hamiltonian must return a ',...
				'Hermitian matrix with the same size as the density matrix.']);
		H1=x;
		n=n+1;
	end
	while n<N
		if strcmpi(varargin{n},'RelTol')
			x=varargin{n+1};
			assert(isrealscalar(x) && x>0,...
				'QuantMech:Solve:LME:InvalidInput',...
				'Input to the relative tolerance must be a positive real scalar.');
			RelTol=x;
		elseif strcmpi(varargin{n},'AbsTol')
			x=varargin{n+1};
			assert(isrealscalar(x) && x>0,...
				'QuantMech:Solve:LME:InvalidInput',...
				'Input to the absolute tolerance must be a positive real scalar.');
			AbsTol=x;
		elseif strcmpi(varargin{n},'InitialStep')
			x=varargin{n+1};
			assert(isrealscalar(x) && x>0 && x<1,...
				'QuantMech:Solve:LME:InvalidInput',...
				'Input to the initial step must be a real scalar between 0 and 1.');
			InitialStep=x;
		elseif strcmpi(varargin{n},'MaxStep')
			x=varargin{n+1};
			assert(isrealscalar(x) && x>0 && x<1,...
				'QuantMech:Solve:LME:InvalidInput',...
				'Input to the maximum step must be a real scalar between 0 and 1.');
			MaxStep=x;
		elseif strcmpi(varargin{n},'Events')
			x=varargin{n+1};
			assert(isa(x,'function_handle') && nargin(x)==2 && nargout(x)==3,...
				'QuantMech:Solve:LME:InvalidInput',...
				['Input to the event function must be a handle to the function ',...
					'which takes two inputs and return three outputs.']);
			Events=x;
		elseif strcmpi(varargin{n},'PlotObs')
			x=varargin{n+1};
			if ~iscell(x)
				x={x};
			end
			assert(~isempty(x) ...
				&& all(cellfun(@(a)iscomplexmatrix(a) && all(size(a)==D),x)),...
				'QuantMech:Solve:LME:InvalidInput',...
				['Input to the output observable must be a complex matrix with ',...
					'the same size as H0 or a cell containing such matrices.']);
			PlotObs=x;
		else
			error('QuantMech:Solve:LME:InvalidInput',...
				'Unrecognized Name-Value pair arguments.');
		end
		n=n+2;
	end
	if n<=N
		warning('QuantMech:Solve:LME:IgnoredInput',...
			'An extra input argument was provided, but is ignored.');
	end
end

%% Solving the LME
x=tspan(end)-tspan(1);
options=odeset(...
	'RelTol',RelTol,...
	'AbsTol',AbsTol,...
	'MaxStep',MaxStep*x);
if InitialStep>0
	options=odeset(options,'InitialStep',InitialStep*x);
end
if ~isempty(PlotObs)
	options=odeset(options,'OutputFcn',@b_plotobs);
end
if isempty(H1)
	options=odeset(options,'Vectorized','on');
end

x=reshape(rho0,D^2,1);
if isempty(Events)
	[t,y]=ode45(@a_odefcn,tspan,[real(x);imag(x)],options);
else
	options=odeset(options,'Events',Events);
	[t,y,te,ye,ie]=ode45(@a_odefcn,tspan,[real(x);imag(x)],options);
	varargout={te,ye,ie};
end

t=t.';
rho=reshape((y(:,1:D^2)+1i*y(:,D^2+1:2*D^2)).',D,D,numel(t));

%% Helper Functions
	function a_dydt=a_odefcn(a_t,a_y)
		a_vrho=a_y(1:D^2,:)+1i*a_y(D^2+1:2*D^2,:);
		a_dvrhodt=L*a_vrho;
		if ~isempty(H1)
			for a_j=1:numel(a_t)
				if H1f
					a_H=H1(a_t(a_j));
				else
					a_H=H1(a_t(a_j),reshape(a_vrho(:,a_j),D,D));
				end
				a_H=a_H/Constant.ReducedPlanck;
				a_dvrhodt(:,a_j)=a_dvrhodt(:,a_j) ...
					+1i*(kron(a_H.',eye(D))-kron(eye(D),a_H))*a_vrho(:,a_j);
			end
		end
		a_dydt=[real(a_dvrhodt);imag(a_dvrhodt)];
	end

	function b_status=b_plotobs(b_t,b_y,b_flag)
		if isempty(b_y)
			b_status=odeplot(b_t,b_y,b_flag);
		else
			b_m=numel(PlotObs);
			b_n=size(b_y,2);
			b_rho=reshape(b_y(1:D^2,:)+1i*b_y(D^2+1:2*D^2,:),D,D,b_n);
			b_evobs=zeros(b_m,b_n);
			for b_i=1:b_m
				for b_j=1:b_n
					b_evobs(b_i,b_j)=abs(trace(b_rho(:,:,b_j)*PlotObs{b_i}));
				end
			end
			b_status=odeplot(b_t,b_evobs,b_flag);
		end
	end
end
