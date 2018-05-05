function s=time2str(t)

if t>259200000
	s=sprintf('%.2f years',t/31536000);
elseif t>60480000
	s=sprintf('%.2f months',t/2592000);
elseif t>8640000
	s=sprintf('%.2f weeks',t/604800);
elseif t>360000
	s=sprintf('%.2f days',t/86400);
elseif t>6000
	s=sprintf('%.2f hours',t/3600);
elseif t>100
	s=sprintf('%.2f minutes',t/60);
else
	s=sprintf('%.2f seconds',t);
end

end