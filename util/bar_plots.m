function [] = bar_plots(features, t0, T, p, K, formats)


[fsz, fszl, lwdt, color, grey, c_true, title_str] = formats{:};

pl = bar(t0:T, features', 1.0, 'stacked', 'FaceColor', 'flat', 'FaceAlpha', 1);
pl(1).CData = color;
pl(2).CData = grey;
hold on
yline(p, 'Color', c_true, 'LineWidth', lwdt)
ylim([0, K])
xlim([t0+1, T])
set(gca, 'FontSize', 15)
legend('Correct', 'Incorrect', 'True Dim', 'FontSize', fszl)
title(title_str, 'FontSize', 20)
ylabel('Number of Features ', 'FontSize', fsz)
xlabel('Time', 'FontSize', fsz)
box on



end