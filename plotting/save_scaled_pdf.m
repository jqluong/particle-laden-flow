function save_scaled_pdf(fname, dims, h)
%Quick script to resize figures. to be edited
set(h, 'Units','inches');
set(h, 'OuterPosition',[4 4 dims(1) dims(2)]);
set(gcf,'PaperPositionMode','auto');

saveTightFigure(h,[fname,'.pdf']);