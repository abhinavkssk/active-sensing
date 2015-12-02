
hist_axes = gca;

h_fig = ancestor(hist_axes,'figure');

% Get x/y limits of axes using axis
limits = axis(hist_axes);




% Cache the original axes position so that axes can be repositioned to
% occupy the space used by the colorstripe if nextplot clears the histogram
% axes.
original_axes_pos = get(hist_axes,'Position');

% In GUIDE, default axes units are characters. In order for axes repositiong
% to behave properly, units need to be normalized.
hist_axes_units_old = get(hist_axes,'units');
set(hist_axes,'Units','Normalized');
% Get axis position and make room for color stripe.
pos = get(hist_axes,'pos');
stripe = 0.075;
set(hist_axes,'pos',[pos(1) pos(2)+stripe*pos(4) pos(3) (1-stripe)*pos(4)])
set(hist_axes,'Units',hist_axes_units_old);

set(hist_axes,'xticklabel','')

% Create axis for stripe
stripe_axes = axes('Parent',get(hist_axes,'Parent'),...
                'Position', [pos(1) pos(2) pos(3) stripe*pos(4)]);
				 				 
limits = axis(stripe_axes);
line(limits([1 2 2 1 1]),limits([3 3 4 4 3]),...
       'LineStyle','-',...
       'Parent',stripe_axes,...
       'Color',get(stripe_axes,'XColor'));

   % Create color stripe
    n=length(x(:,1));
    binInterval = 1/n;
    xdata = [binInterval/2 1-(binInterval/2)];
    %limits(1:2) = range;

        C = (1:n)/n;

    
    % image(X,Y,C) where C is the RGB color you specify. 
    image(xdata,[0 1],repmat(C, [1 1 3]),'Parent',stripe_axes);


%set(stripe_axes,'yticklabel','')
%axis(stripe_axes,limits);

% Put a border around the stripe.



