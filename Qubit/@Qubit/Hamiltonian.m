function H=Hamiltonian(obj)
%% Qubit Hamiltonian
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 17/06/2017
% Last modified: 17/06/2017

H=obj.Energy/2*Pauli.Z;

end
