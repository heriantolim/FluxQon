classdef Optical < Photon
%% Optical Photon
%
% Requires package:
%  - MatCommon_v1.0.0+
%  - PhysConst_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: Microwave.
%
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 10/06/2017
% Last modified: 10/06/2017

methods
	function obj=Optical(varargin)
		obj=obj@Photon(varargin{:});
	end
end

methods (Static=true)
	function x=FrequencyRange()
		x=2*pi*[1.8e14,3.3e14];
	end
end

end
