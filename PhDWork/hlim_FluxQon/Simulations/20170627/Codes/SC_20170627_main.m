% MATLAB R2018a
addpackage('MatCommon_v1.0.0',...
	'MatGraphics_v1.0.0',...
	'PhysConst_v1.0.0',...
	'FluxQon_v1.0.0');

QB=.5+(-3:.01:3)/1e4;
runfunction({'compute','ion','qb','spectrum'},QB);% takes 1 day
runfunction({'compute','qb','mw','spectrum'},QB);% takes 1 day
runfunction({'compute','qon','mw','spectrum'},QB,.91:.03:1.09);% takes 3 weeks

runfunction({'plot','qon','mw','spectrum',1});
runfunction({'plot','qon','mw','spectrum',2});
