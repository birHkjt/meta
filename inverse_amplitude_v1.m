load("C:\Users\86187\Desktop\正向\旧复现\li\library\rod5_2_amp.mat");
%% ============================================================
% 第五层分析：从 matr_idx 反查实际超透镜复振幅
%
% E_rod5_2     : 30 × 100 × 100，相位库
% E_rod5_2_amp : 30 × 100 × 100，振幅库
% matr_idx      : 透镜位置对应的 Rx/Ry 索引
%                 0 表示该位置没有超原子
%
% E_lens = 拼接场
% ============================================================

[Nx, Ny, ~] = size(matr_idx);
Nlambda = size(E_rod5_2, 1);

% 初始化
phase_lens = zeros(Nlambda, Nx, Ny);
amp_lens   = zeros(Nlambda, Nx, Ny);
E_lens     = zeros(Nlambda, Nx, Ny);

% 有效透镜区域
valid_mask = (matr_idx(:,:,1) ~= 0);

% ------------------------------------------------------------
% 从单元库反查每一个透镜位置所使用的超原子
% ------------------------------------------------------------

for ix = 1:Nx
    for iy = 1:Ny

        if ~valid_mask(ix,iy)
            continue;
        end

        irx = matr_idx(ix,iy,1);
        iry = matr_idx(ix,iy,2);

        % 30个波长上的相位
        phase_lens(:,ix,iy) = ...
            E_rod5_2(:,irx,iry);

        % 30个波长上的振幅
        amp_lens(:,ix,iy) = ...
            E_rod5_2_amp(:,irx,iry);

        % 实际复振幅，即拼接场
        E_lens(:,ix,iy) = ...
            amp_lens(:,ix,iy) .* ...
            exp(1i*phase_lens(:,ix,iy));

    end
end

disp('实际超透镜复 pupil（拼接场）重建完成。');


%% ============================================================
% 查看实际选中单元的相位和振幅
% ============================================================

ilambda = 1;     % 第一个波长

figure;
imagesc(squeeze(phase_lens(ilambda,:,:)));
axis image;
colorbar;
title(['Actual phase, \lambda = ',num2str(lambda(ilambda)),' \mum']);

figure;
imagesc(squeeze(amp_lens(ilambda,:,:)));
axis image;
colorbar;
title(['Actual amplitude, \lambda = ',num2str(lambda(ilambda)),' \mum']);


%% ============================================================
% 构造理想 pupil
%
% 理想 pupil 的绝对振幅本身没有物理意义，因此这里先定义为：
%
% E_ideal0 = exp(i*phi_ideal)
%
% 后面计算 overlap 时允许它乘以任意复常数 C。
% ============================================================

R2 = x_grid.^2 + y_grid.^2;

E_ideal = zeros(Nlambda,Nx,Ny);

for ilambda = 1:Nlambda

    % 理想双曲相位
    phi_ideal = ...
        -2*pi/lambda(ilambda) .* ...
        (sqrt(R2 + fl^2) -  sqrt(R^2+fl^2));

    % 只保留实际透镜口径
    phi_ideal(~valid_mask) = 0;

    E_ideal(ilambda,:,:) = ...
        valid_mask .* exp(1i*phi_ideal);

end


%% ============================================================
% Pupil overlap
%
% 这里的标准化 overlap：
%
% eta =
% |<E_ideal,E_lens>|^2 /
% (||E_ideal||^2 ||E_lens||^2)
%
% 因此：
% 1. 不需要人为归一化振幅；
% 2. 整体振幅比例不影响结果；
% 3. 整体相位偏移也不影响结果。
% ============================================================

pupil_overlap = zeros(Nlambda,1);

% 同时计算只考虑相位的 overlap
phase_overlap = zeros(Nlambda,1);

% 相位误差 RMS
phase_rms = zeros(Nlambda,1);

% 相对振幅不均匀程度
amp_std = zeros(Nlambda,1);


for ilambda = 1:Nlambda

    % --------------------------------------------------------
    % 取出二维场
    % --------------------------------------------------------

    E1 = squeeze(E_ideal(ilambda,:,:));
    E2 = squeeze(E_lens(ilambda,:,:));

    phi_actual = squeeze(phase_lens(ilambda,:,:));
    amp_actual = squeeze(amp_lens(ilambda,:,:));

    % --------------------------------------------------------
    % 只在有效透镜区域比较
    % --------------------------------------------------------

    E1 = E1(valid_mask);
    E2 = E2(valid_mask);

    phi_actual = phi_actual(valid_mask);
    amp_actual = amp_actual(valid_mask);

    % ========================================================
    % 1. 完整复振幅 pupil overlap
    % ========================================================

    numerator = abs(sum(conj(E1).*E2))^2;

    denominator = ...
        sum(abs(E1).^2) * ...
        sum(abs(E2).^2);

    pupil_overlap(ilambda) = ...
        numerator / denominator;


    % ========================================================
    % 2. 只考虑相位的 overlap
    % ========================================================

    E_phase_actual = exp(1i*phi_actual);

    numerator_phase = ...
        abs(sum(conj(E1).*E_phase_actual))^2;

    denominator_phase = ...
        sum(abs(E1).^2) * ...
        sum(abs(E_phase_actual).^2);

    phase_overlap(ilambda) = ...
        numerator_phase / denominator_phase;


    % ========================================================
    % 3. 相位误差 RMS
    %
    % 使用 angle(exp(i*delta_phi)) 自动处理 2pi 周期
    % ========================================================

    phi_ideal_valid = angle(E1);

    delta_phi = angle( ...
        exp(1i*(phi_actual - phi_ideal_valid)) );

    phase_rms(ilambda) = ...
        sqrt(mean(delta_phi.^2));


    % ========================================================
    % 4. 相对振幅不均匀程度
    %
    % 由于振幅绝对尺度没有意义，
    % 先除以该波长下实际 pupil 的平均振幅。
    %
    % amp_normalized = amp / mean(amp)
    %
    % 理想值为 1。
    % ========================================================

    amp_normalized = ...
        amp_actual / mean(amp_actual);

    amp_std(ilambda) = ...
        sqrt(mean((amp_normalized - 1).^2));

end


%% ============================================================
% 绘制完整复场 overlap
% ============================================================

figure;

plot(lambda,pupil_overlap,'o-','LineWidth',1.5);

xlabel('\lambda (\mum)');
ylabel('Pupil overlap');
title('Complex pupil overlap');

grid on;


%% ============================================================
% 绘制 phase-only overlap
% ============================================================

figure;

plot(lambda,phase_overlap,'o-','LineWidth',1.5);

xlabel('\lambda (\mum)');
ylabel('Phase-only overlap');
title('Phase-only pupil overlap');

grid on;


%% ============================================================
% 绘制相位 RMS
% ============================================================

figure;

plot(lambda,phase_rms*180/pi,'o-','LineWidth',1.5);

xlabel('\lambda (\mum)');
ylabel('Phase RMS error (deg)');
title('RMS phase error');

grid on;


%% ============================================================
% 绘制相对振幅不均匀程度
% ============================================================

figure;

plot(lambda,amp_std,'o-','LineWidth',1.5);

xlabel('\lambda (\mum)');
ylabel('Normalized amplitude STD');
title('Relative amplitude nonuniformity');

grid on;


%% ============================================================
% 输出结果
% ============================================================

disp('========================================');
disp('第五层分析结果');
disp('========================================');

fprintf('Lambda(um)    Complex       Phase-only      Phase RMS(deg)    Amp STD\n');

for ilambda = 1:Nlambda

    fprintf('%8.4f    %10.6f    %10.6f    %12.4f    %10.6f\n', ...
        lambda(ilambda), ...
        pupil_overlap(ilambda), ...
        phase_overlap(ilambda), ...
        phase_rms(ilambda)*180/pi, ...
        amp_std(ilambda));

end





%% ============================================================
% 不同波长下 y = 0 截线的实际相位与理想相位
% ============================================================

% 选择要绘制的波长
lambda_plot = [5.0, 4.5, 4.0, 3.5];
%lambda_plot = [3.5, 4.0, 4.5, 5];

% 找到 y = 0 对应的网格索引
[~, iy0] = min(abs(y_grid(:,1)));

figure;
hold on;

for ilambda_plot = 1:length(lambda_plot)

    % 找到最接近目标波长的实际数据
    [~, ilambda] = min(abs(lambda - lambda_plot(ilambda_plot)));

    % --------------------------------------------------------
    % 实际拼接场相位
    % --------------------------------------------------------
    E_actual = squeeze(E_lens(ilambda,:,:));

    phi_actual = angle(E_actual);

    % --------------------------------------------------------
    % 理想双曲相位
    % 注意这里使用你修正后的参考：
    %
    % sqrt(R2 + fl^2) - sqrt(R^2 + fl^2)
    % --------------------------------------------------------
    phi_ideal = ...
        -2*pi/lambda(ilambda) .* ...
        (sqrt(R2 + fl^2) - sqrt(R^2 + fl^2));

    % --------------------------------------------------------
    % 取 y = 0 截线
    % --------------------------------------------------------
    valid_line = valid_mask(iy0,:);

x_line = x_grid(iy0,valid_line);

actual_wrapped = phi_actual(iy0,valid_line);
ideal_line = phi_ideal(iy0,valid_line);

% 将实际单元相位放到最接近理想相位的 2pi 分支
actual_line = actual_wrapped + ...
    2*pi*round((ideal_line - actual_wrapped)/(2*pi));

plot(x_line, actual_line, ...
    'LineWidth',1.5, ...
    'DisplayName', ...
    [num2str(lambda(ilambda),'%.3f'),' \mum actual']);

plot(x_line, ideal_line, '--', ...
    'LineWidth',1.2, ...
    'DisplayName', ...
    [num2str(lambda(ilambda),'%.3f'),' \mum ideal']);
end

xlabel('x (m数)');
ylabel('Phase (rad)');
title('Phase profile along y = 0');
legend('Location','best');
grid on;
box on;

hold off;