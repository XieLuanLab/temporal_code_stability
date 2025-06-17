function redrawAcrossDayUnits(frame,currentComposition,shankStr)
im = imread([shankStr 'videoFol/unit' num2str(currentComposition(frame),'%04i') '.png']);
imshow(im)
text(1600,1200,sprintf('%03i',frame),'FontSize',14,'HorizontalAlignment','center','Color','r')
end