classdef Microwave < Photon
%% Microwave Photon
%
% Requires package:
%  - MatCommon_v1.0.0+
%  - PhysConst_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: Optical.
%
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 10/06/2017
% Last modified: 10/06/2017

methods
	function obj=Microwave(varargin)
		obj=obj@Photon(varargin{:});
	end
end

methods (Static=true)
	function x=FrequencyRange()
		x=2*pi*[3e8,3e11];
	end
end

end
