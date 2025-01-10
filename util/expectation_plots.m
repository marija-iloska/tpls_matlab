function [] = expectation_plots(E, time_plot, D,  formats)


[fsz, lwd, col, title_str, y_str, x_str] = formats{:};


for j = 1:D
    plot(time_plot, E(:,j), 'Linewidth', lwd, 'Color', col{j});
    hold on
end
yline(0, 'k', 'linewidth',1.5, 'LineStyle','--')
set(gca, 'FontSize', fsz)
title(title_str, 'FontSize', fsz)
ylabel(y_str, 'FontSize', fsz)
xlabel(x_str, 'FontSize', fsz)
xlim([time_plot(1), time_plot(end)])
box on



end