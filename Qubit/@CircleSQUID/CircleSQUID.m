classdef CircleSQUID < Circle & SQUID
%% Circle SQUID Qubit
%
% Requires package:
%  - MatCommon_v1.0.0+
%  - PhysConst_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: Circle, SQUID.
%
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 15/06/2017
% Last modified: 15/06/2017

methods
	function obj=CircleSQUID(varargin)
		obj=obj@SQUID(varargin{:});
	end
end

end
