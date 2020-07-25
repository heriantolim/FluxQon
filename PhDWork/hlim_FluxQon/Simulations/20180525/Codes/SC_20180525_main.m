% MATLAB R2018a
addpackage('MatCommon_v1.0.0',...
	'MatGraphics_v1.0.0',...
	'PhysConst_v1.0.0',...
	'FluxQon_v1.0.0');

for N=1:3
	for n=1:min(2,N)
		runfunction({'compute','ion','fidelity','vs','g','ratio'},...
			sqrt(N/n)*10.^(-1:.01:1),0,n,N);
		runfunction({'compute','fidelity','vs','detuning'},...
			-10:.1:10,0,n,N,sqrt(N/n));
	end
end% takes a few months

N=2;
K={1:3,4};
CR=[[.5,1,2]*sqrt(2),1];
for n=1:2
	for k=K{n}
		runfunction({'compute','ion','dynamics'},0,0,n,N,1,CR(k));
		runfunction({'compute','qon','dynamics'},0,0,n,N,1,0,CR(k));
	end
end

runfunction({'plot','ion','dynamics',2});
runfunction({'plot','qon','dynamics',2});
runfunction({'plot','efficiency','vs','detuning',1});
runfunction({'plot','efficiency','vs','g','ratio',1});
