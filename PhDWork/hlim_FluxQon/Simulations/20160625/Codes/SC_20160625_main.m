% MATLAB R2015b
% MATLAB R2017a
addpackage('MatCommon_v1.0.0',...
	'MatGraphics_v1.0.0',...
	'PhysConst_v1.0.0',...
	'FluxQon_v1.0.0');

runfunction({'compute','energy','and','state','at','pi'});
runfunction({'compute','energy','near','pi'},pi+(-10:.1:10)/1e4*2*pi);
runfunction({'plot','state','and','energy',1},2/pi);
