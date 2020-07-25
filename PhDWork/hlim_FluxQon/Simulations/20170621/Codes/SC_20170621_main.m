% MATLAB R2018a
addpackage('MatCommon_v1.0.0',...
	'MatGraphics_v1.0.0',...
	'PhysConst_v1.0.0',...
	'FluxQon_v1.0.0');

A=0:179;
runfunction({'compute','energy'},A,0,0);
runfunction({'compute','energy'},A,270,0);
runfunction({'compute','G'},A,0,0);% takes 1 year
runfunction({'compute','G'},A,270,0);% takes 1 year
runfunction({'compute','Grms'},A,0,0,.8,1.6);
runfunction({'compute','Grms'},A,270,0,.8,1.6);

runfunction({'plot','G','distribution',1},1,152,270,0,.8,1.6);
runfunction({'plot','G','histogram',1},152,270,0,.8,1.6);
runfunction({'plot','Grms',1},.8,1.6);

Bmax=runfunction({'compute','Bparmax'});
Brms=runfunction({'compute','Bparrms'});
fprintf('max(Bpar) = %g\n',Bmax(1));
fprintf('max(Bperp) = %g\n',Bmax(2));
fprintf('rms(Bpar) = %g\n',Brms(1));
fprintf('rms(Bperp) = %g\n',Brms(2));
