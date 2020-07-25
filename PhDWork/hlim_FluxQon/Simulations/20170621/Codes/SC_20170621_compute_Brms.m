function B=SC_20170621_compute_Brms()

qb=runfunction({'S','C',20180501,'default','qb'});
ion=runfunction({'S','C',20180501,'default','ion'});

B1=integral2(@(r,z)(Circle.GCFr(qb.Radius,r,z)).^2.*r,...
	0,ion.Radius,0,ion.Height);
B1=sqrt(B1/ion.Height)/ion.Radius;

B2=integral2(@(r,z)(Circle.GCFz(qb.Radius,r,z)).^2.*r,...
	0,ion.Radius,0,ion.Height);
B2=sqrt(2*B2/ion.Height)/ion.Radius;

B=Constant.ReducedPlanck*Constant.VacuumPermeability ...
	*qb.TunnelingFrequency/4/Constant.FluxQuantum*[B1,B2];

end