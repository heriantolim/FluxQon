function p=boltzpdf(E,T)
%% Boltzmann Distribution
%  p=boltzpdf(E,T) returns the probability of a state of energy E being occupied
%    at temperature T.
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
% First created: 17/06/2017
% Last modified: 17/06/2017

assert(isrealarray(E) && all(E(:)>=0),...
	'FluxQon:Distribution:boltzpdf:InvalidInput',...
	'Input to the energy must be a positive real array.');
assert(isrealarray(T) && all(T(:)>=0),...
	'FluxQon:Distribution:boltzpdf:InvalidInput',...
	'Input to the temperature must be a positive real array.');

p=exp(-E./T./Constant.Boltzmann);
p(E~=0 & T==0)=0;
p(E==0 & T==0)=1;

end