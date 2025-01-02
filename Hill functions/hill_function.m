function y = hill_function(x, Vmax, Kd, n)
    % x: 配体浓度
    % Vmax: 最大反应速率
    % Kd: 解离常数
    % n: Hill系数

    % Hill函数计算
    y = (Vmax * x.^n) ./ (Kd^n + x.^n);
end