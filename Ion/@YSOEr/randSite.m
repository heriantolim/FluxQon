function x=randSite(varargin)
%% Random Site
%  x=YSOEr.randSite() returns a random value for the site of an YSOEr, which
%    is either 1 or 2. The probability for each value is 0.5.
%
%  x=YSOEr.randSite(N) or x=YSOEr.randSite(N,false) returns N values of random
%    site in a row vector.
%
%  x=YSOEr.randSite(N,true) returns N values of random site where the
%    distribution of values is balanced, i.e. there is approximately 50% values
%    of 1 and 50% values of 2.
%
% Requires package:
%  - MatCommon_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: randSite.
%
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 06/12/2015
% Last modified: 25/09/2016

%% Input Validation and Parsing
balancing=false;
if nargin==0
	N=1;
else
	N=varargin{1};
	assert(isintegerscalar(N) && N>0,...
		'FluxQon:Ion:YSOEr:randSite:InvalidInput',...
		'Input to the number of elements must be a positive integer scalar.');
	if nargin>1
		balancing=varargin{2};
		assert(isbooleanscalar(balancing),...
			'FluxQon:Ion:YSOEr:randSite:InvalidInput',...
			'Input to the balancing option must be a boolean scalar.');
		if nargin>2
			error('FluxQon:Ion:YSOEr:randSite:WrongNargin',...
				'Incorrect number of input arguments.');
		end
	end
end

%% Main
p=.5;
if balancing
	n2=floor(p*N)+(rand<=p);
	x=[ones(1,n2),2*ones(1,N-n2)];
	x=x(randperm(N));
else
	x=(rand([1,N])>p)+1;
end

end