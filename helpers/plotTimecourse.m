function plotTimecourse(results)
  colors = palettablecolors(size(results,2));
  for i = 1:size(results,2)
    plot(results(:,i),'Color',colors(i,:));
    hold on;
    makepalettable;
  end
  xlabel('time step')
  ylabel('average p(cooperation)');
  hz = 0.5;
  plot(get(gca,'xlim'), [hz hz], 'Color', [0.5,0.5,0.5]);
end