function x=randIsotope(varargin)
%% Random Isotope
%  x=Er.randIsotope() returns a random isotope of Er according to its natural
%    abundance: 22.95% odd (1) and 77.05% even (2).
%
%  x=Er.randIsotope(N) or x=Er.randIsotope(N,false) returns N values of random
%    isotope in a row vector.
%
%  x=Er.randIsotope(N,true) returns N values of random isotope where the
%    distribution of odd and even values is balanced, i.e. there is
%    approximately 22.95% values of 1 and 77.05% values of 2.
%
% Requires package:
%  - Common_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 06/12/2015
% Last modified: 25/09/2016

balancing=false;
if nargin==0
	N=1;
else
	N=varargin{1};
	assert(isintegerscalar(N) && N>0,...
		'FluxQon:Ion:Er:randIsotope:InvalidInput',...
		'Input to the number of elements must be a positive integer scalar.');
	if nargin>1
		balancing=varargin{2};
		assert(isbooleanscalar(balancing),...
			'FluxQon:Ion:Er:randIsotope:InvalidInput',...
			'Input to the balancing option must be a boolean scalar.');
		if nargin>2
			error('FluxQon:Ion:Er:randIsotope:WrongNargin',...
				'Incorrect number of input arguments.');
		end
	end
end

p=.2295;
if balancing
	n2=floor(p*N)+(rand<=p);
	x=[ones(1,n2),2*ones(1,N-n2)];
	x=x(randperm(N));
else
	x=(rand([1,N])>p)+1;
end

end