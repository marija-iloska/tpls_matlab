function [log_p] = log_mvnpdf(y, H, theta, var_y)

    log_p = -0.5/var_y * sum((y - H*theta).^2);

end
