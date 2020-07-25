% MATLAB R2015b
% MATLAB R2017a
% MATLAB R2018a
addpackage('MatCommon_v1.0.0',...
	'MatGraphics_v1.0.0',...
	'PhysConst_v1.0.0',...
	'FluxQon_v1.0.0');

runfunction({'compute','energy','and','state','at','pi'});
runfunction({'compute','frequency','vs','bias'},...
	pi+[-10:.1:-5.1,-5:.02:-3.02,-3:.01:3,3.02:.02:5,5.1:.1:10]/1e4*2*pi);% takes 1 day
runfunction({'compute','frequency','vs','tuning'},...
	[0:.0005:.2,.21:.01:1.79,1.8:.0005:2]*pi);% takes 1 day
runfunction({'plot','potential',1},100);
runfunction({'plot','state','and','frequency',1});