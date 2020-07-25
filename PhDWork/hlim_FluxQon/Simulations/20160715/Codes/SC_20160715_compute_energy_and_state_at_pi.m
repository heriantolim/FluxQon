function SC_20160715_compute_energy_and_state_at_pi(varargin)

P=inputParser;
P.addOptional('JosephsonEnergy',5.203156259013710e-22,@isrealscalar);
P.addOptional('AlphaCoefficient',.7,@isrealscalar);
P.parse(varargin{:});
P=P.Results;

M=100;
qb=a3JJ(P.JosephsonEnergy/10,P.JosephsonEnergy,P.AlphaCoefficient);

[E,W]=qb.solveTISE(M);
savedata(E/P.JosephsonEnergy,'R','energy','at','pi','Rewrite','yes');
savedata(W,'R','state','at','pi','Rewrite','yes');
fprintf('Excitation Frequency = %.6f GHz\n',qb.Frequency/2/pi/1e9);
fprintf('Tunneling Frequency  = %.6f GHz\n',qb.TunnelingFrequency/2/pi/1e9);
fprintf('Frequency Ratio      = %.6f\n',qb.TunnelingFrequency/qb.Frequency);

end