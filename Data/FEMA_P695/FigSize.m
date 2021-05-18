function [  ] = FigSize(W, H )
%% Set figure size
% Size Width Then Height
w = W;
h = H;

set(gcf,'units','inches')

pos = get(gcf,'pos');

set(gcf,'units','inches')
set(gcf,'pos',[pos(1) pos(2) w h ])
set(gcf,'Color','w')

set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [w h]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 w h]);
% set(gcf, 'Position',[0 0 1 1]);

end

