function varargout=LME(tSpan,rho0,L,varargin)
%% Solve the Lindblad Master Equation
%  ...=Solve.LME(tSpan,rho0,L) solves the evolution of a quantum system with a
%  constant Lindbladian superoperator L, given the initial density matrix rho0.
%
%  ...=Solve.LME(...,'Observable',O) computes the values of the observables O
%  throughout the evolution.
%
%  The integration method used to find the solution is the Runge-Kutta(4,5)
%  implemented in the MATLAB ODE Suite. The details can be found in Ref.[1].
%
% Outputs:
%  This function returns one to three variables depending on the number of
%  outputs requested and whether the parameter Observable is specified.
%
%  ========================================
%  | # of outputs | Observable |  Outputs |
%  |   requested  |  specified | returned |
%  ----------------------------------------
%  |            1 | yes        | V        |
%  |            1 | no         | Y        |
%  |            2 | yes        | T, V     |
%  |            2 | no         | T, Y     |
%  |            3 | yes        | T, Y, V  |
%  |            3 | no         | T, Y, [] |
%  ========================================
%
%  T: The times of the solution, returned as a real vector of length N. If the
%     input to tSpan is a vector of length 3 or greater, then T is equal to
%     tSpan. Otherwise, T is automatically generated with bounds set by tSpan.
%
%  Y: The solution density matrix, returned as a complex array of size
%     D-by-D-by-N where D is the length of rho0 and N is the length of T.
%
%  V: The values of the observables, returned only if the Observable parameter
%     is specified. V is either an N-column matrix or a cell array of such
%     matrices. The cell form is given if Observable is a cell array.
%
%  []: An empty numeric array.
%
%  In the output cases where Y is not returned, this function does not keep the
%  solution density matrix in the memory -- considerably improving the
%  computation performance if D is a large number.
%
% Inputs:
%  tSpan: The times at which the solutions are to be found. Specified as a
%         vector of real numbers. If tSpan has two elements or fewer, then the
%         times are automatically generated. For tSpan with two elements, the
%         generated times will start from tSpan(1) and end at tSpan(2) --
%         tSpan(2) must be greater than tSpan(1). For tSpan with one element,
%         the generated times will start from 0 to tSpan(1) -- tSpan(1) must be
%         greater than 0.
%
%  rho0: The initial density matrix, specified as a square complex matrix. It is
%        expected that rho0 has all the properties of a density matrix, i.e.
%        Hermitian, positive semi-definite, and of trace one.
%
%  L: The Lindblad superoperator L, specfied as a complex matrix of size
%     D^2-by-D^2, where D is the size of rho0.
%
% Name-Value Inputs:
%  Additional arguments can be passed to this function in Name-Value syntax.
%
%  Observable: Functions to transform the solution density matrix, specified as
%              a function handle. The function referred by the handle must take
%              an argument that is a complex matrix of size D-by-D and return a
%              numeric vector. This parameter can also be specified as a complex
%              matrix of size D-by-D, or a cell array of function handles and/or
%              D-by-D matrices. When specified as a cell array, the output V is
%              also a cell array where each V{i} is a matrix mapped using the
%              Observable, O{i}. When O{i} is a matrix, the corresponding output
%              is the expectation value of O{i}: V{i}=trace(rho*O{i}).
%
%  RelTol: Relative error tolerance, specified as a positive real scalar.
%          Defaults to 1e-3. At each iteration step, the error relative to the
%          magnitude of the solution is checked if it is below RelTol. If it is
%          not, then the step is repeated with smaller time increments.
%
%  AbsTol: Absolute error tolerance, specified as a positive real scalar.
%          Defaults to 1e-6. The absolute error of the solution is checked at
%          each interation step if it is below AbsTol. If it is not, then the
%          step is repeated with smaller time increments.
%
%  MaxStep: Maximum step size as a fraction of tEnd-tStart, specified as a real
%           scalar between 0 and 1. Defaults to 1e-2. This option sets the upper
%           bound on the size of the time increments.
%
%  DispProgress: Whether to display the computation progress. Specified as
%                either: true, false, 1, 0, 'on', 'off', 'yes', or 'no'.
%                Defaults to false.
%
% Requires package:
%  - MatCommon_v1.0.0+
%
% Tested on:
%  - MATLAB R2017a
%  - MATLAB R2018a
%
% See also: UTE.
%
% References:
%  [1] The MATLAB ODE Suite, L. F. Shampine and M. W. Reichelt, SIAM Journal on
%      Scientific Computing, 18-1, 1997.
%
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 22/04/2018
% Last modified: 26/04/2018

%% Default Parameters
outputV=false;
RelTol=1e-3;
AbsTol=1e-6;
MaxStep=1e-2;
DispProgress=false;

%% Input Parsing
assert(isrealvector(tSpan),...
	'FluxQon:Solve:LME:InvalidInput',...
	'Input to the time interval must be a real vector.');
N=numel(tSpan);
if N>2
	fixedT=true;
	[t0,ix]=min(tSpan);
	tN=max(tSpan);
	columnT=iscolumn(tSpan);
	if columnT
		T=tSpan.';
	else
		T=tSpan;
	end
	tSpan=tN-t0;
else
	if N<2
		assert(tSpan>0,...
			'FluxQon:Solve:LME:InvalidInput',...
			'The time length must be greater than zero.');
		t0=0;
		tN=tSpan;
		columnT=false;
	else
		assert(tSpan(2)>tSpan(1),...
			'FluxQon:Solve:LME:InvalidInput',...
			'The end time must be greater than the start time.');
		t0=tSpan(1);
		tN=tSpan(2);
		columnT=iscolumn(tSpan);
		tSpan=tN-t0;
	end
	fixedT=false;
	ix=1;
	T=t0;
end

D=size(rho0,1);
assert(iscomplexmatrix(rho0) && size(rho0,2)==D,...
	'FluxQon:Solve:LME:InvalidInput',...
	'Input to the initial density matrix must be a square complex matrix.');
if ~ishermitian(rho0)
	warning('FluxQon:Solve:LME:NonidealInput',...
		'Input to the initial density matrix is not Hermitian');
end
if trace(rho0)~=1
	warning('FluxQon:Solve:LME:NonidealInput',...
		'Input to the initial density matrix is not of trace one.');
end
D2=D^2;

assert(iscomplexmatrix(L) && all(size(L)==D2),...
	'FluxQon:Solve:LME:InvalidInput',...
	['Input to the Lindblad superoperator must be a complex matrix with ',...
		'length equal to the squared length of the density matrix.']);

k=1;
K=numel(varargin);
while k<K
	y=varargin{k};
	y1=varargin{k+1};
	if strcmpi(y,'Observable')
		if iscell(y1)
			cellV=true;
			O=y1;
		else
			cellV=false;
			O={y1};
		end
		f1=@(x)isa(x,'function_handle')||(iscomplexmatrix(x)&&all(size(x)==D));
		assert(all(cellfun(f1,O)),...
			'FluxQon:Solve:LME:InvalidInput',...
			['Input to the observables must be either a function handle, ',...
				'a complex matrix of the same size as the density matrix, or ',...
				'a cell containing such function handles and/or matrices.']);
		outputV=true;
	elseif strcmpi(y,'RelTol')
		assert(isrealscalar(y1) && y1>0,...
			'FluxQon:Solve:LME:InvalidInput',...
			'Input to the relative tolerance must be a positive real scalar.');
		RelTol=y1;
	elseif strcmpi(y,'AbsTol')
		assert(isrealscalar(y1) && y1>0,...
			'FluxQon:Solve:LME:InvalidInput',...
			'Input to the absolute tolerance must be a positive real scalar.');
		AbsTol=y1;
	elseif strcmpi(y,'MaxStep')
		assert(isrealscalar(y1) && y1>0 && y1<=1,...
			'FluxQon:Solve:LME:InvalidInput',...
			'Input to the max time step must be a real scalar between 0 and 1.');
		MaxStep=y1;
	elseif strcmpi(y,'DispProgress')
		ME=MException('FluxQon:Solve:LME:InvalidInput',...
			'Invalid input to set the DisplayProgress option.');
		if isstringscalar(y1)
			if any(strcmpi(y1,{'on','off','yes','no'}))
				DispProgress=any(strcmpi(y1,{'on','yes'}));
			else
				throw(ME);
			end
		elseif isscalar(y1)
			if isequal(y1,true)
				DispProgress=true;
			elseif isequal(y1,false)
				DispProgress=false;
			else
				throw(ME);
			end
		else
			throw(ME);
		end
	end
	k=k+2;
end
if k<=K
	warning('FluxQon:Solve:LME:IgnoredInput',...
		'An extra input argument was provided, but is ignored.');
end

%% Output Parsing
if nargout==0
	varargout={};
	return
elseif nargout>3
	error('FluxQon:Solve:LME:WrongNargout',...
		'At most three output arguments can be returned.');
end
outputT=nargout>1;
outputY=~outputV || nargout>2;

%% Initialization
% Display progress.
if DispProgress
	tR='N/A';
	fprintf('Progress: %2d%%. Estimated time remaining: %s.',0,tR);
	tic;
end

% The Runge-Kutta(4,5) method parameters.
pow=1/5;

b11=    1/5;
b21=    3/40;    b22=     9/40;
b31=   44/45;    b32=   -56/15;   b33=   32/9;
b41=19372/6561;  b42=-25360/2187; b43=64448/6561; b44=-212/729;
b51= 9017/3168;  b52=  -355/33;   b53=46732/5247; b54=  49/176; b55=-5103/18656;
b61=   35/384;                    b63=  500/1113; b64= 125/192; b65=-2187/6784;
b66=   11/84;

e1=71/57600; e3=-71/16695; e4=71/1920; e5=-17253/339200; e6=22/525; e7=-1/40;

i12=-183/64;   i13=   37/12;  i14= -145/128;
i32=1500/371;  i33=-1000/159; i34= 1000/371;
i42=-125/32;   i43=  125/12;  i44= -375/64;
i52=9477/3392; i53= -729/106; i54=25515/6784;
i62= -11/7;    i63=   11/3;   i64=  -55/28;
i72=   3/2;    i73=   -4;     i74=    5/2;

% The output variables.
if outputY
	if fixedT
		Y=nan(D,D,N);
	end
	Y(:,:,ix)=rho0;
end
if outputV
	V=cell(size(O));
	K=numel(O);
	typeO=zeros(1,K);
	for k=1:K
		if isnumeric(O{k})
			if fixedT
				V{k}=nan(1,N);
			end
			V{k}(ix)=abs(trace(rho0*O{k}));
		else
			if cellV
				s=sprintf('#%d',k);
			else
				s='\b';
			end
			ME=MException('FluxQon:Solve:LME:InvalidInput',...
				'Invalid function referred to by the observable %s.',s);
			try
				y=O{k}(rho0);
			catch ME1
				ME=addCause(ME,ME1);
				throw(ME);
			end
			if iscolumn(y)
				typeO(k)=1;
			elseif isrow(y)
				typeO(k)=2;
			else
				throw(ME);
			end
			if fixedT
				V{k}=nan(numel(y),N);
			end
			V{k}(:,ix)=y;
		end
	end
end

% The loop variables.
hmin=16*eps(tSpan);
hmax=MaxStep*(tSpan);
t1=t0;
y1=reshape(rho0,D2,1);
f1=L*y1;
n1=norm(y1);
nT=AbsTol/RelTol;
h=max(hmax,.8*RelTol*max(n1,nT)/norm(f1));

%% Integration
% The main loop.
i=1;
while i>0
	% Loop for advancing one step.
	ii=1;
	while ii>0
		% Evaluate the derivatives at times between t1 and t1+h.
		f2=L*(y1+h*(f1*b11));
		f3=L*(y1+h*(f1*b21+f2*b22));
		f4=L*(y1+h*(f1*b31+f2*b32+f3*b33));
		f5=L*(y1+h*(f1*b41+f2*b42+f3*b43+f4*b44));
		f6=L*(y1+h*(f1*b51+f2*b52+f3*b53+f4*b54+f5*b55));
		
		% The estimated solution at time = t1+h.
		t7=t1+h;
		y7=y1+h*(f1*b61+f3*b63+f4*b64+f5*b65+f6*b66);
		f7=L*y7;
		n7=norm(y7);
		
		% Check if the estimated error is within the tolerance.
		e7=h*norm(f1*e1+f3*e3+f4*e4+f5*e5+f6*e6+f7*e7)/max([n1,n7,nT]);
		if e7>RelTol
			% Reject the solution.
			if h>hmin
				% Repeat with a smaller time step.
				if ii<2
					h=max(hmin,h*max(.1,.8*(RelTol/e7)^pow));
				else
					h=max(hmin,.5*h);
				end
				ii=ii+1;
			else
				% Quit the current loop, and then the main loop.
				ME=MException('FluxQon:Solve:LME:FailedToConverge',...
					'The derivatives diverge at time t > %g.',t1);
				ii=-ii;
				i=-i;
			end
		else
			% Accept the solution. Quit the current loop.
			break
		end
	end
	
	if i>0
		if fixedT
			% Interpolate the solution.
			ix=find(T>t1 & T<=t7);
			S=numel(ix);
			if S>0
				s=(T(ix)-t1)/h;
				s2=s.*s;
				y=repmat(y1,1,S)+h*(...
					+f1*(s2.*(i12+s.*(i13+s*i14))+s) ...
					+f3*(s2.*(i32+s.*(i33+s*i34))) ...
					+f4*(s2.*(i42+s.*(i43+s*i44))) ...
					+f5*(s2.*(i52+s.*(i53+s*i54))) ...
					+f6*(s2.*(i62+s.*(i63+s*i64))) ...
					+f7*(s2.*(i72+s.*(i73+s*i74))));
			end
		else
			ix=i+1;
			S=1;
			y=y7;
			if outputT
				T(ix)=t7;
			end
		end
		
		if S>0
			y=reshape(y,D,D,S);
			if outputY
				Y(:,:,ix)=y;
			end
			if outputV
				% Transform the solution.
				k=0;
				while k<K && i>0
					k=k+1;
					if typeO(k)==0
						for s=1:S
							V{k}(ix(s))=abs(trace(y(:,:,s)*O{k}));
						end
					else
						try
							if typeO(k)==1
								for s=1:S
									V{k}(:,ix(s))=O{k}(y(:,:,s));
								end
							else
								for s=1:S
									V{k}(:,ix(s))=O{k}(y(:,:,s)).';
								end
							end
						catch ME1
							if cellV
								s=sprintf('#%d',k);
							else
								s='\b';
							end
							ME=MException('FluxQon:Solve:LME:RuntimeError',...
								['Something went wrong when evaluating ',...
									'the function referred to by the observable %s ',...
									'at time t > %g.'],s,t1);
							ME=addCause(ME,ME1);
							i=-i;
						end
					end
				end
			end
		end
	end

	if i>0
		if t7<tN
			% Adjust the time step, if it has not been.
			if ii==1
				h=min([tN-t7,hmax,h*max(5,.8*(RelTol/e7)^pow)]);
			end

			% Advance the integration one step.
			i=i+1; f1=f7; n1=n7; t1=t7; y1=y7;
			
			% Display progress.
			if DispProgress
				fprintf(repmat('\b',1,32+numel(tR)));
				tR=time2str((tSpan-t7+t0)/(t7-t0)*toc);
				fprintf('%2d%%. Estimated time remaining: %s.',...
					floor((t7-t0)/tSpan*100),tR);
			end
		else
			% Quit the main loop.
			break
		end
	end
end

%% Finalization
if columnT
	if outputT
		T=T.';
	end
	if outputV
		for k=1:K
			V{k}=V{k}.';
		end
	end
end

if outputV && ~cellV
	V=V{1};
end

if outputT
	if outputV
		if nargout<3
			varargout={T,V};
		else
			varargout={T,Y,V};
		end
	else
		if nargout<3
			varargout={T,Y};
		else
			varargout={T,Y,[]};
		end
	end
elseif outputV
	varargout={V};
else
	varargout={Y};
end

if DispProgress
	fprintf(repmat('\b',1,42+numel(tR)));
	tR=time2str(toc);
	if i>0
		fprintf('Complete. ');
	else
		fprintf('Error encountered. ');
	end
	fprintf('Time elapsed: %s.\n',tR);
end

if i<0
	disp(getReport(ME));
end

end