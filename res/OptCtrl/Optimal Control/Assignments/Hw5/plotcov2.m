% PLOTCOV2 - Plots a covariance ellipsoid with axes for a bivariate
%            Gaussian distribution.
%
% Usage:
%   [h, s] = plotcov2(mu, Sigma[, OPTIONS]);
% 
% Inputs:
%   mu    - a 2 x 1 vector giving the mean of the distribution.
%   Sigma - a 2 x 2 symmetric positive semi-definite matrix giving
%           the covariance of the distribution (or the zero matrix).
%
% Options:
%   'conf'      - a scalar between 0 and 1 giving the confidence
%                 interval (i.e., the fraction of probability mass to
%                 be enclosed by the ellipse); default is 0.9.
%   'num-pts'   - if the value supplied is n, then (n + 1)^2 points
%                 to be used to plot the ellipse; default is 20.
%   'label'     - if non-empty, a string that will label the
%                 ellipsoid (default: [])
%   'plot-axes' - a 0/1 flag indicating if the ellipsoid's axes
%                 should be plotted (default: 1)
%   'plot-opts' - a cell vector of arguments to be handed to PLOT3
%                 to contol the appearance of the axes, e.g., 
%                 {'Color', 'g', 'LineWidth', 1}; the default is {}
%   'fill-color' - a color specifier; is this is not [], the
%                  covariance ellipse is filled with this color
%                  (default: [])
% 
% Outputs:
%   h     - a vector of handles on the axis lines
%
% See also: PLOTCOV3

% Copyright (C) 2002 Mark A. Paskin
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
% USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [h, s] = plotcov2(mu, Sigma, varargin)

h = [];
s = [];

if size(Sigma) ~= [2 2], error('Sigma must be a 2 by 2 matrix'); end
if length(mu) ~= 2, error('mu must be a 2 by 1 vector'); end

%Sigma = checkpsd(Sigma);


[p, ...
 n, ...
 label, ...
 plot_axes, ...
 plot_opts, ...
 fill_color] = process_options(varargin, 'conf', 0.9, ...
			       'num-pts', 20, ...
			       'label', [], ...
			       'plot-axes', 1, ...
			       'plot-opts', {}, ...
			       'fill-color', []);
holding = ishold;
% Compute the Mahalanobis radius of the ellipsoid that encloses
% the desired probability mass.
k = conf2mahal(p, 2);
% Scale the covariance matrix so the confidence region has unit
% Mahalanobis distance.
Sigma = Sigma * k;
% The axes of the covariance ellipse are given by the eigenvectors of
% the covariance matrix.  Their lengths (for the ellipse with unit
% Mahalanobis radius) are given by the square roots of the
% corresponding eigenvalues.
[V, D] = eig(full(Sigma));
V = real(V);
D = real(D);
D = abs(D);

% Compute the points on the boundary of the ellipsoid.
t = linspace(0, 2*pi, n);
u = [cos(t(:))'; sin(t(:))'];
w = (V * sqrt(D)) * u;
z = repmat(mu(:), [1 n]) + w;
h = [h; plot(z(1, :), z(2, :), plot_opts{:})];
if (~isempty(fill_color))
  s = patch(z(1, :), z(2, :), fill_color);
end

% Plot the axes.
if (plot_axes)
  hold on;
  L = sqrt(diag(D));
  h = plot([mu(1); mu(1) + L(1) * V(1, 1)], ...
	   [mu(2); mu(2) + L(1) * V(2, 1)], plot_opts{:});
  h = [h; plot([mu(1); mu(1) + L(2) * V(1, 2)], ...
	       [mu(2); mu(2) + L(2) * V(2, 2)], plot_opts{:})];
end


if (~isempty(label))
  th = text(mu(1), mu(2), label);
  set(th, 'FontSize', 18);
  set(th, 'FontName', 'Times');
  set(th, 'FontWeight', 'bold');
  set(th, 'FontAngle', 'italic');
  set(th, 'HorizontalAlignment', 'center');
end

if (~holding & plot_axes) hold off; end

function s = checkpsd(s)

if (any(isnan(s) | isinf(s) | ~isreal(s)))
  warning('S contains complex numbers, Inf, or NaN'); 
end
% Drop any negative eigenvalues.
[V, D] = eig(full(s));
d = real(diag(D));
if (any(d < 0))
  warning(sprintf(['S is not positive semidefinite (min. eig. =' ...
		   ' %0.5g); projecting.'], min(d)));
  d(find(d < 0)) = 0;
  D = diag(d);
  s = V * D * V';
end


% PROCESS_OPTIONS - Processes options passed to a Matlab function.
%                   This function provides a simple means of
%                   parsing attribute-value options.  Each option is
%                   named by a unique string and is given a default
%                   value.
%
% Usage:  [var1, var2, ..., varn[, unused]] = ...
%           process_options(args, ...
%                           str1, def1, str2, def2, ..., strn, defn)
%
% Arguments:   
%            args            - a cell array of input arguments, such
%                              as that provided by VARARGIN.  Its contents
%                              should alternate between strings and
%                              values.
%            str1, ..., strn - Strings that are associated with a 
%                              particular variable
%            def1, ..., defn - Default values returned if no option
%                              is supplied
%
% Returns:
%            var1, ..., varn - values to be assigned to variables
%            unused          - an optional cell array of those 
%                              string-value pairs that were unused;
%                              if this is not supplied, then a
%                              warning will be issued for each
%                              option in args that lacked a match.
%
% Examples:
%
% Suppose we wish to define a Matlab function 'func' that has
% required parameters x and y, and optional arguments 'u' and 'v'.
% With the definition
%
%   function y = func(x, y, varargin)
%
%     [u, v] = process_options(varargin, 'u', 0, 'v', 1);
%
% calling func(0, 1, 'v', 2) will assign 0 to x, 1 to y, 0 to u, and 2
% to v.  The parameter names are insensitive to case; calling 
% func(0, 1, 'V', 2) has the same effect.  The function call
% 
%   func(0, 1, 'u', 5, 'z', 2);
%
% will result in u having the value 5 and v having value 1, but
% will issue a warning that the 'z' option has not been used.  On
% the other hand, if func is defined as
%
%   function y = func(x, y, varargin)
%
%     [u, v, unused_args] = process_options(varargin, 'u', 0, 'v', 1);
%
% then the call func(0, 1, 'u', 5, 'z', 2) will yield no warning,
% and unused_args will have the value {'z', 2}.  This behaviour is
% useful for functions with options that invoke other functions
% with options; all options can be passed to the outer function and
% its unprocessed arguments can be passed to the inner function.

% Copyright (C) 2002 Mark A. Paskin
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
% USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [varargout] = process_options(args, varargin)

% Check the number of input arguments
n = length(varargin);
if (mod(n, 2))
  error('Each option must be a string/value pair.');
end

% Check the number of supplied output arguments
if (nargout < (n / 2))
  error('Insufficient number of output arguments given');
elseif (nargout == (n / 2))
  warn = 1;
  nout = n / 2;
else
  warn = 0;
  nout = n / 2 + 1;
end

% Set outputs to be defaults
varargout = cell(1, nout);
for i=2:2:n
  varargout{i/2} = varargin{i};
end

% Now process all arguments
nunused = 0;
for i=1:2:length(args)
  found = 0;
  for j=1:2:n
    if strcmpi(args{i}, varargin{j})
      varargout{(j + 1)/2} = args{i + 1};
      found = 1;
      break;
    end
  end
  if (~found)
    if (warn)
      warning(sprintf('Option ''%s'' not used.', args{i}));
    else
      nunused = nunused + 1;
      unused{2 * nunused - 1} = args{i};
      unused{2 * nunused} = args{i + 1};
    end
  end
end

% Assign the unused arguments
if (~warn)
  if (nunused)
    varargout{nout} = unused;
  else
    varargout{nout} = cell(0);
  end
end

% CONF2MAHAL - Translates a confidence interval to a Mahalanobis
%              distance.  Consider a multivariate Gaussian
%              distribution of the form
%
%   p(x) = 1/sqrt((2 * pi)^d * det(C)) * exp((-1/2) * MD(x, m, inv(C)))
%
%              where MD(x, m, P) is the Mahalanobis distance from x
%              to m under P:
%
%                 MD(x, m, P) = (x - m) * P * (x - m)'
%
%              A particular Mahalanobis distance k identifies an
%              ellipsoid centered at the mean of the distribution.
%              The confidence interval associated with this ellipsoid
%              is the probability mass enclosed by it.  Similarly,
%              a particular confidence interval uniquely determines
%              an ellipsoid with a fixed Mahalanobis distance.
%
%              If X is an d dimensional Gaussian-distributed vector,
%              then the Mahalanobis distance of X is distributed
%              according to the Chi-squared distribution with d
%              degrees of freedom.  Thus, the Mahalanobis distance is
%              determined by evaluating the inverse cumulative
%              distribution function of the chi squared distribution
%              up to the confidence value.
%
% Usage:
% 
%   m = conf2mahal(c, d);
%
% Inputs:
%
%   c    - the confidence interval
%   d    - the number of dimensions of the Gaussian distribution
%
% Outputs:
%
%   m    - the Mahalanobis radius of the ellipsoid enclosing the
%          fraction c of the distribution's probability mass
%
% See also: MAHAL2CONF

% Copyright (C) 2002 Mark A. Paskin
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
% USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = conf2mahal(c, d)

m = chi2inv(c, d);