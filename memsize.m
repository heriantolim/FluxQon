function memsize(varargin)
%% Estimate the Minimum Memory Size
%  memsize(d) displays the memory size of the Hamiltonian and the Lindblad
%  superoperator given the Hilbert dimension d. d can be specified as an integer
%  vector, in which case the Hilbert dimension is taken as the product of the
%  elements in d.
%
%  memsize(obj1,obj2,...) obtains the Hilbert dimension(s) from the property of
%  obj1, obj2, ...
%
% Requires package:
%  - MatCommon_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: Interaction, Lindblad.
% 
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 19/06/2017
% Last modified: 19/06/2017

if nargin==0
	return
end

if isintegervector(varargin{1})
	M=prod(varargin{1})^2;
else
	M=prod(Hilbert.dimension(varargin{:}))^2;
end
M(2)=M(1)^2;
M=M*16;

fprintf('Memory size for:\n');
fprintf('  Hamiltonian  : %g bytes\n',M(1));
fprintf('  Lindblad     : %g bytes\n',M(2));

end