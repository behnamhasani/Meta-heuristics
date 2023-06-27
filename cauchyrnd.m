
function r= cauchyrnd(varargin)
% USAGE:       r= cauchyrnd(a, b, n, ...)
% 
% Generate random numbers from the Cauchy distribution, r= a + b*tan(pi*(rand(n)-0.5)).
% 
% ARGUMENTS:
% a (default value: 0.0) must be scalars or size(x).
% b (b>0, default value: 1.0) must be scalars or size(x).
% n and onwards (default value: 1) specifies the dimension of the output.
% 
% EXAMPLE:
% r= cauchyrnd(0, 1, 10); % A 10 by 10 array of random values, Cauchy distributed.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Default values
	a=	0.0;
	b=	1.0;
	n=	1;
	
	
	% Check the arguments
	if(nargin >= 1)
		a=	varargin{1};
		if(nargin >= 2)
			b=			varargin{2};
			b(b <= 0)=	NaN;	% Make NaN of out of range values.
			if(nargin >= 3),	n=	[varargin{3:end}];		end
		end
	end
	
	
	% Generate
	r=	cauchyinv(rand(n), a, b);
end
