%处理一般无环带情况的版本
%% 预处理 - 修正版本
%参考文献：18年虞南方消色差
clear; close all; clc;

load("C:\Users\86187\Desktop\正向\旧复现\li\library\rod5_2.mat")
E2=E_rod5_2;
p = 1.65; % 单位微米
m = 68; % 直径上的单元数
R = m * p / 2;
fl = 90;
NA = sin(R/sqrt(R^2+fl^2));
c = 3e8 * 1e6; % μm/s

% 波长和频率定义
lambda = transpose(linspace(5, 3.5, 30)); % 5到3.5μm
frequency = c ./ lambda; % 正确的频率定义
omega = 2*pi*frequency; % 角频率

Rx = transpose(linspace(0.1, 1.6, 100));
Ry = transpose(linspace(0.1, 1.6, 100));

% 相位展开
E2_unwrap = zeros(30, 100, 100); % 此矩阵有正有负
for i = 1:100
    for j = 1:100
        E2_unwrap(:, i, j) = unwrap(E2(:, i, j));
    end
end

% 计算群延迟（对角频率的导数）
group_delay = zeros(30, 100, 100);
for rx = 1:100
    for ry = 1:100
        phase_response = squeeze(E2_unwrap(:, rx, ry));
        % 对角度omega求导
        dphi_domega = gradient(phase_response, omega);
        group_delay(:, rx, ry) = dphi_domega;
    end
end

% 计算平均群延迟（假设线性）
avg_group_delay = squeeze(mean(group_delay, 1)); %正矩阵 100x100

%% 需求计算
[x_grid, y_grid] = meshgrid(linspace(-m/2+0.5, m/2-0.5, m));

% 1. 相位需求（最长波长）
lambda_max = max(lambda);
ideal_phase = 2*pi/lambda_max * (-sqrt(x_grid.^2 + y_grid.^2 + fl^2) + sqrt(R^2+fl^2));%开口向下。
%ideal_phase = 2*pi/lambda_max * (-(x_grid.^2 + y_grid.^2)/(2*fl) + R^2/(2*fl));%二次相位
% 加常数使全为正。即已依照18年虞文章修正

% 2. 群延迟需求（物理意义：光程差/c）
ideal_gd = (-sqrt(x_grid.^2 + y_grid.^2 + fl^2) + sqrt(R^2+fl^2)) / c; %需求为正。应正中最大边缘为0。单位：秒
%ideal_gd = (-(x_grid.^2 + y_grid.^2)/(2*fl) + R^2/(2*fl)) / c;%二次相位

 %imagesc(ideal_gd);
 %%%%

ideal_gd_fs = ideal_gd * 1e15; % 转换为飞秒

%% 库数据准备
% 提取最长波长的相位
phase_at_max_lambda = squeeze(E2_unwrap(1, :, :));

% 库中的群延迟范围
lib_gd_range = [min(avg_group_delay(:)), max(avg_group_delay(:))];
ideal_gd_range = [min(ideal_gd_fs(:)), max(ideal_gd_fs(:))];

fprintf('库群延迟范围: %.2f ~ %.2f fs\n', lib_gd_range(1)*1e15, lib_gd_range(2)*1e15);
fprintf('需求群延迟范围: %.2f ~ %.2f fs\n', 0, ideal_gd_range(2));

%% 改进的筛选算法
matr = zeros(m, m, 2);
matr_idx = zeros(m, m, 2);%matr伴生的索引矩阵。用于数据后处理。
match_type = zeros(m, m); % 记录匹配类型：1=完美匹配，2=良好匹配，3=仅相位匹配，0=透镜外

phase_tol = 2*pi/30; % 相位容差
gd_tol_factor = 0.01; % 群延迟容差系数

count_perfect = 0; % 完美匹配计数
count_good = 0;    % 良好匹配计数
count_phase_only = 0; % 仅相位匹配计数

for i1 = 1:m
    for i2 = 1:m
        x = -m/2*p + (i1-0.5)*p;
        y = -m/2*p + (i2-0.5)*p;
        
        if (x^2 + y^2 > R^2)
            match_type(i1, i2) = 0; % 透镜外
            continue; % 透镜外区域
        end
        
        % 当前点的需求
        target_phase = ideal_phase(i1, i2);
        target_gd = ideal_gd_fs(i1, i2); % 飞秒
        
        % 方法1：优先查找完美匹配
        phase_tol_current = phase_tol;
        gd_tol_current = target_gd * gd_tol_factor;
        
        % 相位筛选（考虑2π包裹）
        phase_diff = mod(phase_at_max_lambda - target_phase + pi, 2*pi) - pi;%范围从(0，2pi)变更为(-pi,pi)
         %其实无所谓了。
        phase_mask = abs(phase_diff) < phase_tol_current;
        
        % 群延迟筛选
        gd_diff = avg_group_delay*1e15 - target_gd; % 转换为飞秒比较
        gd_mask = abs(gd_diff) < gd_tol_current;
        
        % 联合筛选
        joint_mask = phase_mask & gd_mask;
        
        if any(joint_mask(:))
            % 完美匹配：同时满足相位和群延迟
            [idx_rx, idx_ry] = find(joint_mask);
            
            % 选择相位最接近的
            phase_diffs = zeros(length(idx_rx), 1);
            for k = 1:length(idx_rx)
                phase_val = phase_at_max_lambda(idx_rx(k), idx_ry(k));
                phase_diffs(k) = abs(mod(phase_val - target_phase + pi, 2*pi) - pi);
            end
            
            [~, best_idx] = min(phase_diffs);
            matr(i1, i2, 1) = Rx(idx_rx(best_idx));
            matr(i1, i2, 2) = Ry(idx_ry(best_idx));
            
            matr_idx(i1, i2, 1) = idx_rx(best_idx);
            matr_idx(i1, i2, 2) = idx_ry(best_idx);
            
            match_type(i1, i2) = 1; % 完美匹配
            count_perfect = count_perfect + 1;
            
        elseif any(phase_mask(:))
            % 仅相位匹配：放宽群延迟要求
            [idx_rx, idx_ry] = find(phase_mask);
            
            % 计算所有候选的群延迟误差
            gd_errors = zeros(length(idx_rx), 1);
            phase_errors = zeros(length(idx_rx), 1);
            gd_errors_pie = zeros(length(idx_rx), 1);%存放非abs数值，尽量使同上下来减小差距
            
            for k = 1:length(idx_rx)
                gd_errors(k) = abs(avg_group_delay(idx_rx(k), idx_ry(k))*1e15 - target_gd);
                gd_errors_pie(k) = avg_group_delay(idx_rx(k), idx_ry(k))*1e15 - target_gd;
                phase_val = phase_at_max_lambda(idx_rx(k), idx_ry(k));
                phase_errors(k) = abs(mod(phase_val - target_phase + pi, 2*pi) - pi);
            end
            
            % 加权评分：相位误差权重更高
            weights = [0.7, 0.3]; % 相位权重, 群延迟权重
            scores = weights(1)*phase_errors/phase_tol + weights(2)*gd_errors/target_gd;
            %for k = 1:length(idx_rx)
            %    if gd_errors_pie(k) < 0
            %       scores(k)=scores(k)*500;%倾向于选群延迟大于需求而不是小于需求的.感觉这句完全不影响结果
            %    end
            %end

            [~, best_idx] = min(scores);
            matr(i1, i2, 1) = Rx(idx_rx(best_idx));
            matr(i1, i2, 2) = Ry(idx_ry(best_idx));

            matr_idx(i1, i2, 1) = idx_rx(best_idx);
            matr_idx(i1, i2, 2) = idx_ry(best_idx);

            match_type(i1, i2) = 2; % 良好匹配
            count_good = count_good + 1;
            
        else
            % 放宽相位要求
            phase_tol_expanded = phase_tol * 3;
            phase_mask = abs(phase_diff) < phase_tol_expanded;
            
            if any(phase_mask(:))
                [idx_rx, idx_ry] = find(phase_mask);
                
                % 选择相位最接近的
                phase_diffs = zeros(length(idx_rx), 1);
                for k = 1:length(idx_rx)
                    phase_val = phase_at_max_lambda(idx_rx(k), idx_ry(k));
                    phase_diffs(k) = abs(mod(phase_val - target_phase + pi, 2*pi) - pi);
                end
                
                [~, best_idx] = min(phase_diffs);
                matr(i1, i2, 1) = Rx(idx_rx(best_idx));
                matr(i1, i2, 2) = Ry(idx_ry(best_idx));

                matr_idx(i1, i2, 1) = idx_rx(best_idx);
                matr_idx(i1, i2, 2) = idx_ry(best_idx);

                match_type(i1, i2) = 3; % 仅相位匹配
                count_phase_only = count_phase_only + 1;
            else
                % 完全无法匹配，选择最接近的相位
                [~, min_idx] = min(abs(phase_diff(:)));
                [idx_rx, idx_ry] = ind2sub(size(phase_diff), min_idx);
                matr(i1, i2, 1) = Rx(idx_rx);
                matr(i1, i2, 2) = Ry(idx_ry);
                
                matr_idx(i1, i2, 1) = idx_rx;
                matr_idx(i1, i2, 2) = idx_ry;
                
                match_type(i1, i2) = 3; % 归为仅相位匹配
            end
        end
    end
end

%% 结果统计
fprintf('\n========== 匹配结果统计 ==========\n');
fprintf('完美匹配（相位+群延迟）: %d/%d (%.1f%%)\n', ...
    count_perfect, m*m, count_perfect/(m*m)*100);
fprintf('良好匹配（相位优先）: %d/%d (%.1f%%)\n', ...
    count_good, m*m, count_good/(m*m)*100);
fprintf('仅相位匹配: %d/%d (%.1f%%)\n', ...
    count_phase_only, m*m, count_phase_only/(m*m)*100);
fprintf('透镜外区域: %d/%d (%.1f%%)\n', ...
    sum(match_type(:)==0), m*m, sum(match_type(:)==0)/(m*m)*100);

%% 可视化验证
figure('Position', [400, 10, 800, 800]);

% 1. 设计的相位分布
subplot(2,2,1);
phase_map = zeros(m,m);
for i1 = 1:m
    for i2 = 1:m
        if matr(i1,i2,1) > 0
            rx_idx = find(Rx == matr(i1,i2,1));
            ry_idx = find(Ry == matr(i1,i2,2));
            phase_map(i1,i2) = phase_at_max_lambda(rx_idx, ry_idx);
        end
    end
end
imagesc(phase_map);
colormap('turbo');
title('设计的相位分布');
colorbar;
axis equal tight;

% 2. 设计的群延迟分布
subplot(2,2,2);
gd_map = zeros(m,m);
for i1 = 1:m
    for i2 = 1:m
        if matr(i1,i2,1) > 0
            rx_idx = find(Rx == matr(i1,i2,1));
            ry_idx = find(Ry == matr(i1,i2,2));
            gd_map(i1,i2) = avg_group_delay(rx_idx, ry_idx) * 1e15;
        end
    end
end
imagesc(gd_map);
title('设计的群延迟分布 (fs)');
colorbar;
axis equal tight;

% 3. 误差分析
subplot(2,2,3);
phase_error = mod(phase_map - ideal_phase + pi, 2*pi) - pi;
imagesc(abs(phase_error));
title('相位误差分布 (rad)');
colorbar;
axis equal tight;

% 4. 匹配类型可视化
subplot(2,2,4);

% 创建自定义颜色映射
match_colors = [
    0.8 0.8 0.8;  % 灰色 - 透镜外 (0)
    0.2 0.8 0.2;  % 绿色 - 完美匹配 (1)
    1.0 0.8 0.2;  % 橙色 - 良好匹配 (2)
    0.8 0.2 0.2;  % 红色 - 仅相位匹配 (3)
];

% 创建匹配类型图
match_type_image = zeros(m, m, 3);
for i = 1:m
    for j = 1:m
        type_idx = match_type(i, j) + 1; % +1因为索引从1开始
        match_type_image(i, j, :) = match_colors(type_idx, :);
    end
end

image(match_type_image);
title('单元匹配类型');
axis equal tight;

% 添加图例
hold on;
% 创建图例手动元素
legend_labels = {'透镜外', '完美匹配', '良好匹配', '仅相位匹配'};
for i = 1:4
    plot(NaN, NaN, 's', 'MarkerFaceColor', match_colors(i, :), ...
        'MarkerEdgeColor', 'k', 'MarkerSize', 10, 'DisplayName', legend_labels{i});
end
legend('Location', 'southoutside', 'Orientation', 'horizontal', 'NumColumns', 2);

%% 中轴线上群延迟对比
figure('Position', [100, 100, 800, 400]);

% 计算中轴线索引（y=0对应的行）
if mod(m, 2) == 0
    center_row = m/2 + 1;  % 偶数时取中心偏上的行
else
    center_row = (m+1)/2;  % 奇数时取中心行
end

% 提取中轴线上的X坐标
x_positions = linspace(-R, R, m);  % X坐标（μm）

% 1. 理论需求曲线
theory_gd_mid = ideal_gd_fs(center_row, :);  % 理论群延迟（飞秒）

% 2. 实际设计曲线（从最终设计中提取）
actual_gd_mid = zeros(1, m);
for i = 1:m
    if matr(center_row, i, 1) > 0  % 有效单元
        rx_idx = find(Rx == matr(center_row, i, 1));
        ry_idx = find(Ry == matr(center_row, i, 2));
        actual_gd_mid(i) = avg_group_delay(rx_idx, ry_idx) * 1e15;  % 转换为飞秒
    else
        actual_gd_mid(i) = NaN;  % 无效单元标记为NaN
    end
end

% 绘制对比
plot(x_positions, theory_gd_mid, 'r-', 'LineWidth', 1, 'DisplayName', '理想需求');
hold on;
plot(x_positions, actual_gd_mid, 'b--o', 'LineWidth', 0.5, 'DisplayName', '实际采样');
hold off;

% 美化图形
grid on;
box on;
xlabel('X 位置 (μm)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('群延迟 (fs)', 'FontSize', 12, 'FontWeight', 'bold');
title('中轴线(y=0)群延迟对比', 'FontSize', 14, 'FontWeight', 'bold');
%legend('Location', 'northwest', 'FontSize', 11);

% 添加透镜边缘标记
hold on;
lens_radius = R;
plot([-lens_radius, -lens_radius], ylim, 'k:', 'LineWidth', 1);
plot([lens_radius, lens_radius], ylim, 'k:', 'LineWidth', 1);
%text(-lens_radius, max(ylim)*0.95, ' 透镜边缘', 'FontSize', 9, 'VerticalAlignment', 'top');
hold off;

% 添加简要统计
%rms_error = sqrt(mean((actual_gd_mid - theory_gd_mid).^2));
%text(0.02, 0.98, sprintf('RMS误差: %.1f fs', rms_error), ...
%    'Units', 'normalized', 'FontSize', 10, 'VerticalAlignment', 'top');
%% 基于下一节的群延迟探讨
omegamax=c*2*pi/3.5;
omegamin=c*2*pi/5;
deltaomega=omegamax-omegamin;

avg_group_delay0=zeros(100,100);
for i1=1:100
    for i2=1:100
deltap=E2_unwrap(30,i1,i2)-E2_unwrap(1,i1,i2);
avg_group_delay0(i1,i2)=deltap/deltaomega;
    end
end
%这个处理之后相比原来的处理会稍微宽一点，不过没有本质区别，可能只是边缘问题。
%% 目标：库的群延迟分析图样
E2_unwrap0=E2_unwrap+pi;
E2_unwrap1=squeeze(E2_unwrap0(1,:,:));%最长波刚好覆盖（0，2Pi）不超过
gpg=zeros(10000,2);
omegamax=c*2*pi/3.5;
omegamin=c*2*pi/5;
deltaomega=omegamax-omegamin;
for i1=1:100
    for i2=1:100
        gpg(100*(i1-1)+i2,1)=E2_unwrap1(i1,i2)/pi;
        gpg(100*(i1-1)+i2,2)=deltaomega*avg_group_delay0(i1,i2)/pi;
    end
end
%% 

figure('Position', [200, 200, 800, 400]);
scatter(gpg(:,1),gpg(:,2),1,"b",".");%库图
xlim([0,2])
ylim([0,4.5])

hold on

%% ===========================================================
% GD-Phase library 与透镜中轴需求轨迹
% ============================================================

figure('Position',[200,200,800,400]);

% ---------- 1. Library GD-Phase 点云 ----------
phase_lib = mod(phase_at_max_lambda, 2*pi) / pi;
gd_lib = deltaomega * avg_group_delay0 / pi;

gpg = [phase_lib(:), gd_lib(:)];

scatter(gpg(:,1), gpg(:,2), 1, 'b', '.');
hold on;

% ---------- 2. 实际单元中心坐标 ----------
x_cell = -R + p/2 : p : R-p/2;

% ---------- 3. 理想需求 ----------
% 使用与主程序完全相同的定义：
% phase：5 um
% GD：光程差 / c

phase_req = 2*pi/lambda_max .* ...
    (-sqrt(x_cell.^2 + fl^2) + sqrt(R^2 + fl^2));

gd_req = ...
    (-sqrt(x_cell.^2 + fl^2) + sqrt(R^2 + fl^2)) / c;

% ---------- 4. 转换到 GD-Phase 图坐标 ----------
phase_req = mod(phase_req,2*pi) / pi;
gd_req = deltaomega * gd_req / pi;

% ---------- 5. 只取一侧，避免左右完全重复 ----------
half_idx = 1:ceil(length(x_cell)/2);

scatter(phase_req(half_idx), ...
        gd_req(half_idx), ...
        10, ...
        'filled');

xlabel('相位 (\pi)');
ylabel('群延迟 (\Delta\omega GD/\pi)');
xlim([0,2]);
ylim([0,4.5]);

grid on;
box on;

hold off;
xlabel("相位（Pi）")
ylabel("群延迟(频差相移)")
hold off


