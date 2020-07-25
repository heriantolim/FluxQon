% MATLAB R2018a
addpackage('MatCommon_v1.0.0',...
	'MatGraphics_v1.0.0',...
	'PhysConst_v1.0.0',...
	'FluxQon_v1.0.0');

N=1e4;
for n=1:2
	runfunction({'compute','ion','dynamics'},0,0,n,N,1,sqrt(N/n));
	runfunction({'compute','qon','dynamics'},0,0,n,N,1,0,sqrt(N/n));
end

runfunction({'plot','ion','dynamics',1});
runfunction({'plot','qon','dynamics',1});
