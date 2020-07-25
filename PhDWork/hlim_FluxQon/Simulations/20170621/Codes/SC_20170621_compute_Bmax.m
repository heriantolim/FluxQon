function B=SC_20170621_compute_Bmax()

N=4001;

qb=runfunction({'S','C',20180501,'default','qb'});
ion=runfunction({'S','C',20180501,'default','ion'});

B1=0;
for r=linspace(0,ion.Radius,N)
	for z=linspace(0,ion.Height,N)
		B1=max(B1,Circle.GCFr(qb.Radius,r,z));
	end
end

B2=2*pi/qb.Radius;

B=Constant.ReducedPlanck*Constant.VacuumPermeability ...
	*qb.TunnelingFrequency/4/Constant.FluxQuantum*[B1,B2];

end