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
% Inputs:
%  - fMax : Maximum probe frequency to be resolved, specified as a positive real
%           scalar.
%  - tMax : Maximum probe time, specified as a positive real scalar.
%  - L    : Lindblad superoperator, specified as a complex square matrix.
%  - R    : Readout operator, specified as a square complex matrix.
%  - d    : Subspace dimensions, specified as a positive integer vector. The
%           squared product of the elements in d must equal to the dimension of
%           L.
%  - n    : Index of the photon subspace, specified as a positive integer
%           scalar.
%
% Outputs:
%  - S    : Radiation spectrum, returned as a two-row matrix. The first row is
%           the frequency (x-data points), and the second row is the intensity
%           (y-data points) of the radiation signals.
%
% Requires package:
%  - Common_v1.0.0+
%  - QuantMech_v1.0.0+
%
% Tested on:
%  - MATLAB R2017a
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 17/06/2017
% Last modified: 17/06/2017

assert(isrealscalar(fMax) && fMax>0,...
	'FluxQon:Readout:RadiationSpectrum:InvalidInput',...
	'Input to the maximum probe frequency must be a positive real scalar.');
assert(isrealscalar(tMax) && tMax>0,...
	'FluxQon:Readout:RadiationSpectrum:InvalidInput',...
	'Input to the maximum probe time must be a positive real scalar.');

D=size(L);
assert(iscomplexmatrix(L) && diff(D)==0,...
	'FluxQon:Readout:RadiationSpectrum:InvalidInput',...
	'Input to the Lindblad superoperator must be a complex square matrix.');
D=sqrt(D(1));
assert(isintegerscalar(D),...
	'FluxQon:Readout:RadiationSpectrum:InvalidInput',...
	'The dimension of the Lindblad superoperator must be an integer squared.');

if nargin>6
	error('FluxQon:Readout:RadiationSpectrum:WrongNargin',...
		'Incorrect number of input arguments.');
elseif nargin==4
	assert(iscomplexmatrix(varargin{1}) && all(size(varargin{1})==D),...
		'FluxQon:Readout:RadiationSpectrum:InvalidInput',...
		['Input to the readout operator in this syntax must be a square ',...
			'complex matrix whose dimension is equal to the square-root of ',...
			'the dimension of the Lindblad superoperator.']);
	R=varargin{1};
else
	assert(isintegervector(varargin{1}) && all(varargin{1}>0),...
		'FluxQon:Readout:RadiationSpectrum:InvalidInput',...
		'Input to the subspace dimensions must be a positive integer vector.');
	assert(prod(varargin{1})==D,...
		'FluxQon:Readout:RadiationSpectrum:InvalidInput',...
		['The squared total dimension of the subspaces must equal to ',...
			'the dimension of the Lindblad superoperator.']);
	d=varargin{1};

	assert(isintegerscalar(varargin{2}) && varargin{2}>0,...
		'FluxQuon:Readout:RadiationSpectrum:InvalidInput',...
		'Input to the subspace index must be a positive integer scalar.');
	n=varargin{2};

	if nargin==5
		R=Annihilation(d(n));
	elseif isa(varargin{3},'Photon')
		R=varargin{3}.Annihilation;
	elseif iscomplexmatrix(varargin{3}) && all(size(varargin{3})==d(n))
		R=varargin{3};
	else
		error('FluxQon:Readout:RadiationSpectrum:InvalidInput',...
			['Input to the readout operator in this syntax must be either ',...
				'a Photon object or a square complex matrix whose dimension ',...
				'is equal to the dimension of the corresponding subspace.']);
	end
	R=Operator.kron(d,n,R);
end

% Steady-state density matrix
fprintf('Radiation Spectrum\n');
fprintf('------------------\n');
fprintf('Finding steady state solution : ');
rho_ss=null(L);
J=size(rho_ss,2);
if J==0
	error('FluxQon:Readout:RadiationSpectrum:NoSolution',...
		'The steady-state solution does not exist.');
elseif J==1
	fprintf('complete.\n');
else
	fprintf('\n');
	warning('FluxQon:Readout:RadiationSpectrum:NonUniqueSolution',...
		'The nullity of the Lindblad superoperator is greater than 1.');
end

% Sampling points
t=0:1/2/fMax:tMax;
J=numel(t);
fprintf('Number of sampling points     : %g.\n',J);

% Correlation function of the radiated photons
S=zeros(1,J);
rho_ss=reshape(R*reshape(rho_ss(:,1),D,D),D^2,1);
R=R';
fprintf('Computing correlation function:  ');
for j=1:J
	fprintf('%2d%%',floor(100*j/J));
	S(j)=trace(R*reshape(expm(L*t(j))*rho_ss,D,D));
	fprintf('\b\b\b');
end
fprintf('\b\bcomplete.\n');

% Radiation spectrum
J=2^nextpow2(J);
f=0:2*fMax/J:fMax;
S=fft(S,J);
S=abs(S(1:J/2+1))/J;
S(2:end-1)=2*S(2:end-1);
S=[f;S];

end
