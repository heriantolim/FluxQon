function memsize(varargin)
%% Estimate the Minimum Memory Size
%  memsize(obj1,obj2,...,objN) prints the memory size of the Hamiltonian and the
%    Lindblad superoperator if constructed from obj1, obj2, ..., objN.
%
% Requires package:
%  - Common_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: Interaction, Lindblad.
% 
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 19/06/2017
% Last modified: 19/06/2017

if nargin==0
	return
end

M=prod(Hilbert.dimension(varargin{:}))^2;
M(2)=M(1)^2;
M=M*16;

fprintf('Memory size for:\n');
fprintf('  Hamiltonian  : %g Bytes\n',M(1));
fprintf('  Lindblad     : %g Bytes\n',M(2));

end