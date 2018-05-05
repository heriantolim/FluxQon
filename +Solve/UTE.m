function varargout=UTE(tSpan,psi0,H0,varargin)
%% Solve the Unitary Time Evolution
%  ...=Solve.UTE(tSpan,psi0,H0) solves the unitary evolution of a quantum system
%  with a constant Hamiltonian H0, given the initial state psi0.
%
%  ...=Solve.UTE(...,'Observable',O) computes the values of the observables O
%  throughout the evolution.
%
%  Solving a unitary evolution with a time-varying Hamiltonian is not supported.
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
%  Y: The solution state, returned as a complex matrix of size D-by-N where D is
%     the length of psi0 and N is the length of T.
%
%  V: The values of the observables, returned only if the Observable parameter
%     is specified. V is either an N-column matrix or a cell array of such
%     matrices. The cell form is given if Observable is a cell array.
%
%  []: An empty numeric array.
%
%  In the output cases where Y is not returned, this function does not keep the
%  solution state in the memory -- considerably improving the computation
%  performance if D is a large number.
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
%  psi0: The initial state, specified as a complex vector.
%
%  H0: The constant Hamiltonian, specfied as a Hermitian matrix of size D-by-D,
%      where D is the length of psi0.
%
% Name-Value Inputs:
%  Additional arguments can be passed to this function in Name-Value syntax.
%
%  Observable: Functions to transform the solution state, specified as a
%              function handle. The function referred by the handle must take an
%              argument that is a complex vector of size D-by-1 and return a
%              numeric vector. This parameter can also be specified as a complex
%              matrix of size D-by-D, or a cell array of function handles and/or
%              D-by-D matrices. When specified as a cell array, the output V is
%              also a cell array where each V{i} is a matrix mapped using the
%              Observable, O{i}. When O{i} is a matrix, the corresponding output
%              is the expectation value of O{i}: V{i}=psi.'*O{i}*psi.
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
% See also: LME.
%
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 22/04/2018
% Last modified: 26/04/2018

%% Default Parameters
outputV=false;
DefaultStep=1e-2;% times the smallest oscillation period.
DispProgress=false;

%% Input Parsing
assert(isrealvector(tSpan),...
	'FluxQon:Solve:UTE:InvalidInput',...
	'Input to the time interval must be a real vector.');
N=numel(tSpan);
if N>2
	[t0,ix]=min(tSpan);
	columnT=iscolumn(tSpan);
	if columnT
		T=tSpan.'-t0;
	else
		T=tSpan-t0;
	end
	tSpan=max(tSpan)-t0;
else
	if N<2
		assert(tSpan>0,...
			'FluxQon:Solve:UTE:InvalidInput',...
			'The time length must be greater than zero.');
		t0=0;
		columnT=false;
	else
		assert(tSpan(2)>tSpan(1),...
			'FluxQon:Solve:UTE:InvalidInput',...
			'The end time must be greater than the start time.');
		t0=tSpan(1);
		columnT=iscolumn(tSpan);
		tSpan=tSpan(2)-t0;
	end
	ix=1;
end

assert(iscomplexvector(psi0),...
	'FluxQon:Solve:UTE:InvalidInput',...
		'Input to the initial state must be a complex vector.');
D=numel(psi0);
n=norm(psi0);
if n~=1
	warning('FluxQon:Solve:UTE:NonidealInput',...
		'The input initial state is not normalized. It will be made so.');
	psi0=psi0/n;
end

assert(iscomplexmatrix(H0) && all(size(H0)==D) && ishermitian(H0),...
	'FluxQon:Solve:UTE:InvalidInput',...
	['Input to the Hamiltonian must be a Hermitian matrix of ',...
		'the same length as the initial state.']);

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
		assert(all(cellfun(@(x)isa(x,'function_handle') ...
				|| (iscomplexmatrix(x) && all(size(x)==D)),O)),...
			'FluxQon:Solve:UTE:InvalidInput',...
			['Input to the observables must be either a function handle, ',...
				'a complex matrix of the same length as the initial state, or ',...
				'a cell containing such function handles and/or matrices.']);
		outputV=true;
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
	warning('FluxQon:Solve:UTE:IgnoredInput',...
		'An extra input argument was provided, but is ignored.');
end

%% Output Parsing
if nargout==0
	varargout={};
	return
elseif nargout>3
	error('FluxQon:Solve:UTE:WrongNargout',...
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
	n=2;
end

% The output variables.
if outputY
	Y=nan(D,N);
	Y(:,ix)=psi0;
end
if outputV
	V=cell(size(O));
	K=numel(O);
	typeO=zeros(1,K);
	for k=1:K
		if isnumeric(O{k})
			V{k}=nan(1,N);
			V{k}(ix)=abs(psi0'*O{k}*psi0);
		else
			if cellV
				s=sprintf('#%d',k);
			else
				s='\b';
			end
			ME=MException('FluxQon:Solve:UTE:InvalidInput',...
				'Invalid function referred to by the observable %s.',s);
			try
				y=O{k}(psi0);
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
			V{k}=nan(numel(y),N);
			V{k}(:,ix)=y;
		end
	end
end

% The loop variables.
[H0,W]=eig(-1i*H0/Constant.ReducedPlanck);% H0 is now the eigenvectors.
W=diag(W).';
if N<3
	N=ceil(max(abs(W))/2/pi/DefaultStep*tSpan);
	T=linspace(0,tSpan,N);
end

%% Evaluation
I=1:N;
I(ix)=[];
ME=[];
for i=I
	y=((H0.*exp(W*T(i)))/H0)*psi0;
	if outputY
		Y(:,i)=y;
	end
	if outputV
		for k=1:K
			if typeO(k)==0
				V{k}(i)=abs(y'*O{k}*y);
			else
				try
					if typeO(k)==1
						V{k}(:,i)=O{k}(y);
					else
						V{k}(:,i)=O{k}(y).';
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
							'at time t = %g.'],s,T(i));
					ME=addCause(ME,ME1);
					break
				end
			end
		end
	end
	if DispProgress && isempty(ME) && n<N
		fprintf(repmat('\b',1,32+numel(tR)));
		tR=time2str((N-n)/n*toc);
		fprintf('%2d%%. Estimated time remaining: %s.',floor(n/N*100),tR);
		n=n+1;
	end
end

%% Finalization
if outputT
	T=T+t0;
end
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
	if isempty(ME)
		fprintf('Complete. ');
	else
		fprintf('Error encountered. ');
	end
	fprintf('Time elapsed: %s.\n',tR);
end

if ~isempty(ME)
	disp(getReport(ME));
end

end