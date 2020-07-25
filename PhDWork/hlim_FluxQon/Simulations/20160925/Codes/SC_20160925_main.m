% MATLAB R2015b
% MATLAB R2017a
% MATLAB R2018a
addpackage('MatCommon_v1.0.0',...
	'MatGraphics_v1.0.0',...
	'PhysConst_v1.0.0',...
	'FluxQon_v1.0.0');

runfunction({'compute','energy','vs','rot1angle'},0:180);
runfunction({'plot','energy','vs','rot1angle',1});
