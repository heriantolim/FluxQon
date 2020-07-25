% MATLAB R2018a
addpackage('MatCommon_v1.0.0',...
	'MatGraphics_v1.0.0',...
	'PhysConst_v1.0.0',...
	'FluxQon_v1.0.0');

for n=[1,3,4]
	runfunction({'compute','qon','mw','dynamics'},0,n,2,.5);
end
for n=1:2
	runfunction({'compute','ion','mw','dynamics'},0,n,2,.06459);
end
for n=3:4
	runfunction({'compute','ion','mw','dynamics'},0,n,2,.06369);
end

n=2;
for dt=[0,.5]
	for cr=[.5,1,2]
		runfunction({'compute','qon','mw','dynamics'},1,n,0,dt,cr);
		runfunction({'compute','qon','mw','dynamics'},1,n,2,dt,cr);
	end
end
for dt=[-.5,-.25,.25,1]
	runfunction({'compute','qon','mw','dynamics'},0,n,0,dt);
end
for dt=[.25,1,2]
	runfunction({'compute','qon','mw','dynamics'},0,n,2,dt);
end

N=1:4;
DT=-.88:.02:5;
runfunction({'compute','fidelity','vs','detuning'},N,DT,0);
runfunction({'compute','fidelity','vs','detuning'},N,DT,2);

n=2;
CR=10.^(-1:.01:1);
for dt=[0,.5,1]
	for cm=0:3
		runfunction({'compute','fidelity','vs','g','ratio'},n,CR,cm,dt);
	end
end

DT=-.25:.001:.25;
runfunction({'compute','osc','spectrum','vs','detuning'},2,DT,'ComponentID',{1});
DT=-.88:.01:5;
runfunction({'compute','osc','spectrum','vs','detuning'},2,DT,'ComponentID',{2});
CT=.0001*(1:500);
runfunction({'compute','osc','spectrum','vs','concentration'},2,CT,'ComponentID',{1});

runfunction({'plot','qon','mw','dynamics',3},2);
runfunction({'plot','qon','mw','dynamics',4},2);
runfunction({'plot','dynamics',1},4);
runfunction({'plot','max','qb','excitation','vs','detuning',1});
runfunction({'plot','exchange','fidelity','vs','detuning',2});
runfunction({'plot','exchange','fidelity','vs','g','ratio',1},0);
runfunction({'plot','exchange','fidelity','comparison',1});
runfunction({'plot','exchange','time','vs','detuning',1});
runfunction({'plot','osc','spectrum','vs','detuning',2},2);
runfunction({'plot','osc','spectrum','vs','concentration',2},2);
runfunction({'plot','g','effective','vs','tuning','flux',1});
