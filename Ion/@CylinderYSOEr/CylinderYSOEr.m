classdef CylinderYSOEr < Cylinder & YSOEr
%% Er3+ in a Cubic Y2SiO5 Host
%
% Requires package:
%  - MatCommon_v1.0.0+
%  - PhysConst_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: Cylinder, YSOEr.
%
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 15/06/2016
% Last modified: 15/06/2017

methods
	function obj=CylinderYSOEr(varargin)
		obj=obj@YSOEr(varargin{:});
	end
end

end
