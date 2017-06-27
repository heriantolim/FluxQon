function H=Hamiltonian(obj)
%% Qubit Hamiltonian
%
% Requires package:
%  - QuantMech_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 17/06/2017
% Last modified: 17/06/2017

H=obj.Energy/2*Pauli.Z;

end