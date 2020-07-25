function SC_20160625_compute_energy_and_state_at_pi(varargin)

P=inputParser;
P.addOptional('JosephsonEnergy',Constant.FluxQuantum/2/pi*1e-6,@isrealscalar);
P.parse(varargin{:});
P=P.Results;

M=1000;
qb=SQUID(P.JosephsonEnergy/25,P.JosephsonEnergy*2/pi,P.JosephsonEnergy);

[E,W]=qb.solveTISE(M);
savedata(E,'R','energy','at','pi','Rewrite','yes');
savedata(W,'R','state','at','pi','Rewrite','yes');
fprintf('Excitation Frequency = %.6f GHz\n',qb.Frequency/2/pi/1e9);
fprintf('Tunneling Frequency  = %.6f GHz\n',qb.TunnelingFrequency/2/pi/1e9);
fprintf('Frequency Ratio      = %.6f\n',qb.TunnelingFrequency/qb.Frequency);

end