clear all
close all
clc

R = 1000;
var_y = 0.001;   % Variance
ps = 4;     % Sparsity percent
dy = 7;      % System dimension


str_dy = num2str(dy);
str_k = num2str(dy - ps);
str_v = num2str(var_y);
str_R = num2str(R);

all_T = [30, 60, 120, 240, 480];


for t = 1:length(all_T)

    % Get correct Data
    T = all_T(t);
    str_T = num2str(T);

    % COnstruct full filename to load
    filename = join(['Results/T', str_T, '_K', str_dy, '_k', str_k, '_v', str_v, ...
        '_R', str_R, '.mat']);

    % Load data
    load(filename)

    % Variables of interest
    time_mcmc = time_mcmc(:,1);
    time_orls = time_orls(:,1);
    run_time_ratio(t) = mean(time_mcmc./time_orls);
    avg_time_mcmc(t) = mean(time_mcmc);
    avg_time_orls(t) = mean(time_orls);

    mcmc{t} = mcmc_run(:,1);
    orls{t} = orls_run(:,1);

    % Number of times chosen correctly 
    mcmc_first(t) = sum(mcmc{t} == 1);
    orls_first(t) = sum(orls{t} == 1);

    mcmc_second(t) = sum(mcmc{t} == 2);
    orls_second(t) = sum(orls{t} == 2);


    figure;
    subplot(2,1,1)
    histogram(orls{t}, 'FaceColor', [176, 123, 173]/256, 'FaceAlpha', 0.8, ...
        'EdgeColor', [80, 0, 110]/256, 'EdgeAlpha', 1, 'LineWidth', 1.5)
    orls_x = get(gca, 'xTick');
    xticks(unique(round(orls_x)));
    ylim([0,R])
    xlim([0,5.5])
    hold on
    set(gca, 'FontSize', 15)
    title('pjORLS', 'FontSize', 15)
    ylabel('Runs', 'FontSize', 15)
    xlabel('Rank of Correct Model', 'FontSize', 15)
    grid on


    subplot(2,1,2)
    histogram(mcmc{t}, 'FaceColor', [9, 173, 168]/256, 'FaceAlpha', 0.8, ...
        'EdgeColor', [31, 61, 60]/256, 'EdgeAlpha', 1, 'LineWidth', 1.5)
    mcmc_x = get(gca, 'xTick');
    xticks(unique(round(mcmc_x)));
    ylim([0,R])
    xlim([0,5.5])
    set(gca, 'FontSize', 15)
    title('rjMCMC', 'FontSize', 15)
    ylabel('Runs', 'FontSize', 15)
    xlabel('Rank of Correct Model', 'FontSize', 15)
    grid on



    figure;
    histogram(orls{t}, 'FaceColor', [176, 123, 173]/256, 'FaceAlpha', 0.5, ...
        'EdgeColor', [80, 0, 110]/256, 'EdgeAlpha', 1)
    orls_x = get(gca, 'xTick');
    xticks(unique(round(orls_x)));
    ylim([0,R])
    hold on
    histogram(mcmc{t}, 'FaceColor', [9, 173, 168]/256, 'FaceAlpha', 0.3, ...
        'EdgeColor', [31, 61, 60]/256, 'EdgeAlpha', 1)
    set(gca, 'FontSize', 15)
    ylabel('Runs', 'FontSize', 15)
    xlabel('Rank of Correct Model', 'FontSize', 15)
    legend('pjORLS', 'rjMCMC', 'FontSize',15)
    grid on





end

figure;
plot(all_T, run_time_ratio)


figure;
plot(all_T, orls_first, 'linewidth', 1)
hold on
plot(all_T, mcmc_first, 'linewidth', 1)
hold on
plot(all_T, orls_second, 'linewidth', 1, 'LineStyle','--')
hold on
plot(all_T, mcmc_second, 'linewidth', 1, 'LineStyle','--')

