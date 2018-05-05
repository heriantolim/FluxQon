function S=RadiationSpectrum(fMax,tMax,L,varargin)
%% Radiation Spectrum
%  S=RadiationSpectrum(fMax,tMax,L,R) computes the radiation spectrum of a
%    quantum system described by the Lindblad superoperator L. The spectrum is
%    measured by the readout operator R, and is probed from frequency 0 to fMax
%    for a time span from 0 to tMax. The squared dimension of R must equal to
%    the dimension of L.
%
%  S=RadiationSpectrum(fMax,tMax,L,d,n) sets the readout operator as an
%    annihilation operator in the subspace at index n. The dimensions of the
%    subspace and the other subspaces are listed in the integer vector d.
%
%  S=RadiationSpectrum(fMax,tMax,L,d,n,R) sets the readout operator to R, which
%    acts only on the subspace at index n. The dimension of R must equal to the
%    dimension of the subspace listed in d.
%
%  S=RadiationSpectrum(fMax,tMax,L,d,n,obj), where obj is a Photon object, sets
%    the readout operator as the annihilation operator of the object.
%
% Outputs:
%  - S: Radiation spectrum, returned as a two-row matrix. The first row is the
%       frequency (x-data points), and the second row is the intensity (y-data
%       points) of the radiation signals.
%
% Inputs:
%  - fMax: The maximum probe frequency, specified as a positive real scalar.
%  - tMax: The maximum probe time, specified as a positive real scalar.
%  - L: The Lindblad superoperator, specified as a complex square matrix.
%  - R: The readout operator, specified as a square complex matrix.
%  - d: The subspace dimensions, specified as a positive integer vector. The
%       squared product of the elements in d must equal to the dimension of L.
%  - n: Index of the photon subspace, specified as a positive integer scalar.
%
% Name-Value Inputs:
%  Additional arguments can be passed to this function in Name-Value syntax.
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
%
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 17/06/2017
% Last modified: 17/06/2017

%% Default Parameters
DispProgress=false;

%% Input Parsing
assert(isrealscalar(fMax) && fMax>0,...
	'FluxQon:Readout:RadiationSpectrum:InvalidInput',...
	'Input to the maximum probe frequency must be a positive real scalar.');
assert(isrealscalar(tMax) && tMax>0,...
	'FluxQon:Readout:RadiationSpectrum:InvalidInput',...
	'Input to the maximum probe time must be a positive real scalar.');

D=size(L);
assert(iscomplexmatrix(L) && D(1)==D(2),...
	'FluxQon:Readout:RadiationSpectrum:InvalidInput',...
	'Input to the Lindblad superoperator must be a complex square matrix.');
D=sqrt(D(1));
assert(isintegerscalar(D),...
	'FluxQon:Readout:RadiationSpectrum:InvalidInput',...
	'The dimension of the Lindblad superoperator must be a squared integer.');

ME=MException('FluxQon:Readout:RadiationSpectrum:WrongNargin',...
	'Incorrect number of input arguments.');
k=1;
K=numel(varargin);
if K==0
	throw(ME);
elseif isintegervector(varargin{k}) && all(varargin{k}>0)
	% Input syntax #1.
	d=varargin{k};
	assert(prod(d)==D,...
		'FluxQon:Readout:RadiationSpectrum:InvalidInput',...
		['The squared total dimension of the subspaces must equal to ',...
			'the dimension of the Lindblad superoperator.']);
	k=k+1;
	
	if K<2
		throw(ME);
	elseif isintegerscalar(varargin{k}) && varargin{k}>0
		n=varargin{k};
	else
		error('FluxQuon:Readout:RadiationSpectrum:InvalidInput',...
			'Input to the subspace index must be a positive integer scalar.');
	end
	k=k+1;
	
	if K==2 || isstringscalar(varargin{k})
		R=Annihilation(d(n));
	else
		if isa(varargin{k},'Photon')
			R=varargin{k}.Annihilation;
		elseif iscomplexmatrix(varargin{k}) && all(size(varargin{k})==d(n))
			R=varargin{k};
		else
			error('FluxQon:Readout:RadiationSpectrum:InvalidInput',...
				['Input to the readout operator in this syntax must be either ',...
					'a Photon object or a square complex matrix whose dimension ',...
					'is equal to the dimension of the corresponding subspace.']);
		end
		k=k+1;
	end
	R=Operator.kron(d,n,R);
elseif iscomplexmatrix(varargin{k}) && all(size(varargin{k})==D)
	% Input syntax #2.
	R=varargin{k};
	k=k+1;
else
	error('FluxQon:Readout:RadiationSpectrum:InvalidInput',...
		['The fourth input argument must be either a positive integer vector ',...
			'or a complex square matrix whose dimension is equal to the ',...
			'square-root of the dimension of the Lindblad superoperator.']);
end

while k<K
	W=varargin{k+1};
	if strcmpi(varargin{k},'DispProgress')
		ME=MException('FluxQon:Solve:LME:InvalidInput',...
			'Invalid input to set the DisplayProgress option.');
		if isstringscalar(W)
			if any(strcmpi(W,{'on','off','yes','no'}))
				DispProgress=any(strcmpi(W,{'on','yes'}));
			else
				throw(ME);
			end
		elseif isscalar(W)
			if isequal(W,true)
				DispProgress=true;
			elseif isequal(W,false)
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

%% Main
if DispProgress
	tR='N/A';
	fprintf('Progress: %2d%%. Estimated time remaining: %s.',0,tR);
	tic;
end

% The steady-state density matrix.
K=null(L);
N=size(K,2);
if N==0
	error('FluxQon:Readout:RadiationSpectrum:NoSolution',...
		'The steady-state solution does not exist.');
elseif N>1
	warning('FluxQon:Readout:RadiationSpectrum:NonUniqueSolution',...
		'The nullity of the Lindblad superoperator is greater than 1.');
end

% The sampling points.
T=0:1/2/fMax:tMax;
N=numel(T);

% Eigendecomposition of L.
[L,W]=eig(L);% L is now the eigenvectors.
W=diag(W).';

% Correlation function of the radiated photons.
S=zeros(1,N);
K=reshape(R*reshape(K(:,1),D,D),D^2,1);
R=R';
for n=1:N
	S(n)=trace(R*reshape(((L.*exp(W*T(n)))/L)*K,D,D));
	if DispProgress && n<N
		fprintf(repmat('\b',1,32+numel(tR)));
		tR=time2str((N-n)/n*toc);
		fprintf('%2d%%. Estimated time remaining: %s.',floor(n/N*100),tR);
	end
end

% Radiation spectrum.
N=2^nextpow2(N);
F=0:2*fMax/N:fMax;
S=fft(S,N);
S=abs(S(1:N/2+1))/N;
S(2:end-1)=2*S(2:end-1);
S=[F;S];

if DispProgress
	fprintf(repmat('\b',1,42+numel(tR)));
	tR=time2str(toc);
	fprintf('Complete. Time elapsed: %s.\n',tR);
end

end
