function H=Hamiltonian(obj)
%% Photon Hamiltonian
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

if obj.PumpEnergy==0 || obj.PumpPower==0
	H=obj.Energy*obj.Number;
else
	H=(obj.Energy-obj.PumpEnergy)*obj.Number ...
		+obj.PumpPower*(obj.Annihilation+obj.Creation);
end

end