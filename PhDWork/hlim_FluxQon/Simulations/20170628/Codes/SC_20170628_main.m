% MATLAB R2018a
addpackage('MatCommon_v1.0.0',...
	'MatGraphics_v1.0.0',...
	'PhysConst_v1.0.0',...
	'FluxQon_v1.0.0');

CouplingRatio=[.5,1,2];
IB=.85:.0005:1.15;
for R=CouplingRatio% takes 9 days
	runfunction({'compute','qb','mw','spectrum'},.5,R);
	runfunction({'compute','qon','mw','spectrum'},IB,.5,R);
end

runfunction({'plot','qon','mw','spectrum',1});
