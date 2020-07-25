% MATLAB R2018a
addpackage('MatCommon_v1.0.0',...
	'MatGraphics_v1.0.0',...
	'PhysConst_v1.0.0',...
	'FluxQon_v1.0.0');

CouplingRatio=[.5,1,2];
D=0:.001:.5;
for R=CouplingRatio% takes 7 days
	runfunction({'compute','qon','mw','spectrum'},D,1,R);
end

runfunction({'plot','qon','mw','spectrum',1});
