function [] = WhitePcshow()

% Function to paint white the pcshow botton
set(gcf,'color','w'); % fonde de la figura a blanco
set(gca,'color','w'); % fondo de la caja de la nube a blanco
set(gca, 'XColor', [1 1 1], 'YColor', [1 1 1], 'ZColor', [1 1 1]); % ejes a blanco

end

