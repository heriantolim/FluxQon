% MATLAB R2018a
addpackage('MatCommon_v1.0.0',...
	'MatGraphics_v1.0.0',...
	'PhysConst_v1.0.0',...
	'FluxQon_v1.0.0');

n=2;
DT=[0,.5];
CR=[.5,1,2];

runfunction({'compute','qb','mw','dynamics'},0,n,2);
runfunction({'compute','ion','mw','dynamics'},0,n,2);
runfunction({'compute','qon','dynamics'},0,n,2,1);
runfunction({'compute','qon','dynamics'},0,n,2,2);

for dt=DT
	for cr=CR
		runfunction({'compute','qon','mw','dynamics'},1,n,0,dt,cr);
	end
end

for dt=DT
	runfunction({'compute','qon','mw','dynamics'},1,n,2,dt);% takes 10 days
end

n=2;
runfunction({'plot','qon','dynamics',1},n);
runfunction({'plot','qon','mw','dynamics',1},n);
runfunction({'plot','qon','mw','dynamics',2},n);
runfunction({'plot','qon','mw','dynamics',3},n);
runfunction({'show','fidelity','stats',1},n);
