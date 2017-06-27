classdef CubeYSOEr < Cube & YSOEr
%% Er3+ in a Cubic Y2SiO5 Host
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
% See also: Cube, YSOEr.
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 15/06/2016
% Last modified: 15/06/2017

methods
	function obj=CubeYSOEr(varargin)
		obj=obj@YSOEr(varargin{:});
	end
end

end
