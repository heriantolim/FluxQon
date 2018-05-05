function f=fidelity(p1,p2)
%% Fidelity
%  f=fidelity(p1,p2) computes the fidelity of the states or density matrices p1
%  and p2. Both p1 and p2 must have the same length. If one of p1 and p2 is a
%  state and the other is a matrix, then the state will be converted to a
%  density matrix via p*p'.
%
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 03/08/2017
% Last modified: 27/04/2018

d1=size(p1);
d2=size(p2);
tf1=any(d1==1);
tf2=any(d2==1);
assert(iscomplexmatrix(p1) && iscomplexmatrix(p2) ...
		&& (tf1 || d1(1)==d1(2)) ...
		&& (tf2 || d2(1)==d2(2)),...
	'FluxQon:State:fidelity:InvalidInput',...
	'The inputs must be a vector or a square matrix of complex numbers.');

d1=max(d1);
d2=max(d2);
assert(d1==d2,...
	'FluxQon:State:fidelity:InvalidInput',...
	'Both inputs must have the same length');

if tf1 && tf2
	% Both inputs are states.
	f=abs(p1'*p2)^2;
	return
end

% Convert states to density matrices.
if tf1
	p1=p1*p1';
end
if tf2
	p2=p2*p2';
end

% Compute the fidelity.
if isdiag(p1)
	p1=sqrt(p1);
else
	p1=sqrtm(p1);
end
p2=p1*p2*p1;
if isdiag(p2)
	p2=sqrt(p2);
else
	p2=sqrtm(p2);
end
f=abs(trace(p2))^2;

end