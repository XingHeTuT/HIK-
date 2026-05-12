clc; clear; close all;

%% ====== 需要你修改的路径 ======
inputFolder  = 'D:\LiLinhan_2026\相位反演\HIK相位测量\更改焦平面步长\Zemax_txt';
outputFolder = 'D:\LiLinhan_2026\相位反演\HIK相位测量\更改焦平面步长\Zemax_xlsx';
% ==============================

if ~exist(inputFolder, 'dir'), error('输入文件夹不存在：%s', inputFolder); end
if ~exist(outputFolder, 'dir'), mkdir(outputFolder); end

files = dir(fullfile(inputFolder, '*.txt')); 
fprintf('找到 %d 个文件。\n', numel(files));

for i = 1:numel(files)
    txtPath = fullfile(inputFolder, files(i).name);
    [~, baseName, ~] = fileparts(files(i).name);
    xlsxPath = fullfile(outputFolder, baseName + ".xlsx");
    
    fprintf('------------------------------------------------------\n');
    fprintf('[%d/%d] 正在读取：%s\n', i, numel(files), files(i).name);
    
    try
        % 1. 读取数据 (使用最稳健的 readmatrix 文本模式)
        M = readmatrix(txtPath, 'FileType', 'text');
        
        [rows_old, cols_old] = size(M);
        fprintf('    原始尺寸：%d x %d\n', rows_old, cols_old);
        
        % ================= 自动清洗逻辑 =================
        if ~isempty(M)
            % --- A. 检测并删除第一列 ---
            % 逻辑：如果不为空，且满足以下任一条件则删除
            % 1. 差分为1 (序号)
            % 2. 几乎全为0 (杂项)
            % 3. 全是NaN (读取错误的文本列)
            if size(M, 2) > 1
                col1 = M(:, 1);
                d = diff(col1);
                
                is_sequence = ~isempty(d) && all(abs(d(~isnan(d)) - 1) < 1e-4);
                is_zeros    = all(abs(col1(~isnan(col1))) < 1e-4);
                is_nans     = all(isnan(col1));
                
                if is_sequence || is_zeros || is_nans
                    M(:, 1) = []; 
                    fprintf('    -> [清洗] 检测到第一列为序号/杂项，已删除。\n');
                end
            end
            
            % --- B. 检测并删除第一行 ---
            if size(M, 1) > 1
                row1 = M(1, :);
                d = diff(row1);
                
                is_sequence = ~isempty(d) && all(abs(d(~isnan(d)) - 1) < 1e-4);
                is_zeros    = all(abs(row1(~isnan(row1))) < 1e-4);
                is_nans     = all(isnan(row1));
                
                if is_sequence || is_zeros || is_nans
                    M(1, :) = [];
                    fprintf('    -> [清洗] 检测到第一行为序号/杂项，已删除。\n');
                end
            end
        end
        % ==============================================
        
        [rows_new, cols_new] = size(M);
        fprintf('    最终内存尺寸：%d x %d', rows_new, cols_new);
        if rows_new == 512 && cols_new == 640
            fprintf(' (符合预期)\n');
        else
            fprintf('\n');
        end
        
        % >>>>> 关键修复：写入前强制删除旧文件 <<<<<
        if exist(xlsxPath, 'file')
            try
                delete(xlsxPath);
                fprintf('    -> [文件操作] 旧 Excel 文件已删除，准备写入新文件。\n');
            catch
                warning('无法删除旧文件，请确保 Excel 已关闭！');
            end
        end
        
        % 保存
        writematrix(M, xlsxPath);
        fprintf('    -> 保存成功：%s\n', xlsxPath);
        
    catch ME
        warning('处理失败：%s\n错误信息：%s', files(i).name, ME.message);
    end
end
fprintf('------------------------------------------------------\n');
fprintf('全部完成。请检查 Excel 文件尺寸。\n');