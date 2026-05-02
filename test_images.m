clear all;
% get data
input_folder_images = 'DIV2K_valid_HR_data';
images = {'0871'};
name_folder = "results_2026-05-02_1314";
folder = 'results/img/' + name_folder;
if ~exist(char(folder), "dir")
    mkdir(char(folder));
end
transforms_file = "results/results_2026-05-02_1314/results_1.mat";
config_file = "results/results_2026-05-02_1314/global_settings.mat";
load(transforms_file);
load(config_file);
n = patch_size_list(1);
T0 = T0_list(1);
blocks = [];
for i = 1:length(images)
    filepath = fullfile(input_folder_images, [images{i}, '.mat']);
    data = struct2cell(load(filepath));
    img = data{1};
    blocks = [blocks, my_im2col(img, [sqrt(n), sqrt(n)], sqrt(n))];
end
br = mean(blocks); 
YH = blocks - (ones(n, 1) * br); 
STY = T0 * ones(1, size(YH, 2));
% get transforms
W0 = kron(dctmtx(sqrt(n)), dctmtx(sqrt(n)));
T_sc = results(1).structured_conditioned.transform * W0;
T_cond = results(1).unstructured_conditioned.transform;
T_bresler = results(1).bresler_method.transform * W0;
T_dct = W0; % Definire matrice DCT

% apply tranform
X1 = T_sc * YH;
[s] = sort(abs(X1), 'descend'); 
X_sc = X1 .* (bsxfun(@ge, abs(X1), s(STY)));

X1 = T_cond * YH;
[s] = sort(abs(X1), 'descend'); 
X_cond = X1 .* (bsxfun(@ge, abs(X1), s(STY)));

X1 = T_bresler * YH;
[s] = sort(abs(X1), 'descend'); 
X_bresler = X1 .* (bsxfun(@ge, abs(X1), s(STY)));

X1 = T_dct * YH;
[s] = sort(abs(X1), 'descend'); 
X_dct = X1 .* (bsxfun(@ge, abs(X1), s(STY)));

% reconstruct the image
T_sc_inv = inv(T_sc);
T_cond_inv = inv(T_cond);
T_bresler_inv = inv(T_bresler);
T_dct_inv = inv(T_dct);

Y_hat_sc = T_sc_inv * X_sc;
Y_hat_sc = Y_hat_sc + (ones(size(Y_hat_sc,1),1) * br);

Y_hat_cond = T_cond_inv * X_cond;
Y_hat_cond = Y_hat_cond + (ones(size(Y_hat_cond,1),1) * br);

Y_hat_bresler = T_bresler_inv * X_bresler;
Y_hat_bresler = Y_hat_bresler + (ones(size(Y_hat_bresler,1),1) * br);

Y_hat_dct = T_dct_inv * X_dct;
Y_hat_dct = Y_hat_dct + (ones(size(Y_hat_dct,1),1) * br);

% save results
imageSize = [1024, 1024];  
patchSize = [sqrt(n), sqrt(n)];     
numPatchesPerImage = (imageSize(1) / patchSize(1)) * (imageSize(2) / patchSize(2));  
startIdx = 1;
for i = 1:length(images)
    endIdx = startIdx + numPatchesPerImage - 1;
    filepath = fullfile(input_folder_images, [images{i}, '.mat']);
    data = struct2cell(load(filepath));
    img = data{1};
    img_db = uint8(img);
    
    % Bresler
    Y_hat_bresler_img = Y_hat_bresler(:, startIdx:endIdx);
    I_hat_bresler = col2im(Y_hat_bresler_img, patchSize, imageSize, 'distinct');
    I_hat_bresler_db = uint8(I_hat_bresler);
    fprintf('Image: %s | Method: Bresler | PSNR: %.2f | MSE: %.2f | SSIM: %.4f\n', ...
        images{i}, psnr(I_hat_bresler_db, img_db), immse(I_hat_bresler_db, img_db), ssim(I_hat_bresler_db, img_db));
    imwrite(uint8(I_hat_bresler), sprintf('results/img/' + name_folder + '/%s_%s.png', "bresler", images{i}));
    
    % Unstructured Conditioned
    Y_hat_cond_img = Y_hat_cond(:, startIdx:endIdx);
    I_hat_cond = col2im(Y_hat_cond_img, patchSize, imageSize, 'distinct');
    I_hat_cond_db = uint8(I_hat_cond);
    fprintf('Image: %s | Method: Unstructured Conditioned | PSNR: %.2f | MSE: %.2f | SSIM: %.4f\n', ...
        images{i}, psnr(I_hat_cond_db, img_db), immse(I_hat_cond_db, img_db), ssim(I_hat_cond_db, img_db));
    imwrite(uint8(I_hat_cond), sprintf('results/img/' + name_folder + '/%s_%s.png', "unstructured_conditoned", images{i}));
    
    % Structured Conditioned
    Y_hat_sc_img = Y_hat_sc(:, startIdx:endIdx);
    I_hat_sc = col2im(Y_hat_sc_img, patchSize, imageSize, 'distinct');
    I_hat_sc_db = uint8(I_hat_sc);
    fprintf('Image: %s | Method: Structured Conditioned | PSNR: %.2f | MSE: %.2f | SSIM: %.4f\n', ...
        images{i}, psnr(I_hat_sc_db, img_db), immse(I_hat_sc_db, img_db), ssim(I_hat_sc_db, img_db));
    imwrite(uint8(I_hat_sc), sprintf('results/img/' + name_folder + '/%s_%s.png', "structured_conditioned", images{i}));

    % DCT
    Y_hat_dct_img = Y_hat_dct(:, startIdx:endIdx);
    I_hat_dct = col2im(Y_hat_dct_img, patchSize, imageSize, 'distinct');
    I_hat_dct_db = uint8(I_hat_dct);
    fprintf('Image: %s | Method: DCT | PSNR: %.2f | MSE: %.2f | SSIM: %.4f\n', ...
        images{i}, psnr(I_hat_dct_db, img_db), immse(I_hat_dct_db, img_db), ssim(I_hat_dct_db, img_db));
    imwrite(uint8(I_hat_dct), sprintf('results/img/' + name_folder + '/%s_%s.png', "dct", images{i}));
    
    startIdx = endIdx + 1;
end