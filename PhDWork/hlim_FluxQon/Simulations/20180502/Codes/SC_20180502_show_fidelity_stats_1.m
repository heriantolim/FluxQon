function [A,B]=SC_20180502_show_fidelity_stats_1(n)

A=[];
B=[];
i=0;
for k=1:3
	S=loaddata('R','plot','qon','mw','dynamics',k,n);
	for l=1:size(S.Y,1)
		i=i+1;
		A(i)=mean(S.Y{l,2}{2});
		B(i)=min(S.Y{l,2}{2});
	end
end

fprintf('Average fidelity = 1 - %g\n',1-mean(A));
fprintf('Minimum fidelity = 1 - %g\n',1-min(B));

end