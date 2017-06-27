classdef Circle3JJ < Circle & a3JJ
%% Circle 3JJ Flux Qubit
%
% Requires package:
%  - Common_v1.0.0+
%  - PhysConst_v1.0.0+
%  - QuantMech_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: Circle, a3JJ.
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 15/06/2017
% Last modified: 15/06/2017

methods
	function obj=Circle3JJ(varargin)
		obj=obj@a3JJ(varargin{:});
	end
end

end
