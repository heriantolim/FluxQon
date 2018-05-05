classdef SQUID < FluxQubit
%% SQUID Qubit
%
% Constructor syntax:
%  obj=SQUID() creates a SQUID object with default properties.
%
%  obj=SQUID(wz,wx) creates a SQUID object with a frequency of wz and a
%    tunneling frequency of wx.
%
%  obj=SQUID(EC,EL,EJ) creates a SQUID object with a Coulomb energy EC, an
%    inductance energy EL, and a Josephson energy EJ. The energy and tunneling
%    energy are computed upon calling obj.solveTISE.
%
%  obj=SQUID('PropertyName',PropertyValue,...) sets the properties of the object
%    in Name-Value pair syntax.
%
% Requires package:
%  - MatCommon_v1.0.0+
%  - PhysConst_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%  - MATLAB R2017a
%
% See also: FluxQubit.
%
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 13/12/2015
% Last modified: 15/06/2017

methods
	function obj=SQUID(varargin)
		obj=obj@FluxQubit(varargin{:});
	end

	[E,W,dX]=solveTISE(obj,N)
end

methods (Access=protected)
	function autoSetEnergy(obj)
		try
			obj.solveTISE();
		catch ME1
			if isempty(regexpi(ME1.identifier,':MissingData$','once'))
				rethrow(ME1);
			end
		end
	end

	function autoSetTunnelingEnergy(obj)
		try
			obj.solveTISE();
		catch ME1
			if isempty(regexpi(ME1.identifier,':MissingData$','once'))
				rethrow(ME1);
			end
		end
	end
end

end
