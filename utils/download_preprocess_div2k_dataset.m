% download data from http://data.vision.ee.ethz.ch/cvl/DIV2K/DIV2K_valid_HR.zip

url = 'http://data.vision.ee.ethz.ch/cvl/DIV2K/DIV2K_valid_HR.zip';
filename = 'DIV2K_valid_HR.zip';

websave(filename, url);
unzip(filename);
delete(filename);

input_folder = 'DIV2K_valid_HR';     
output_folder_data = 'DIV2K_valid_HR_data_full';
% output_folder_png = 'DIV2K_valid_HR_png';
target_size = [1024, 1024];

if ~exist(output_folder_data, 'dir')
    mkdir(output_folder_data);
end
% if ~exist(output_folder_png, 'dir')
%     mkdir(output_folder_png);
% end

files = dir(fullfile(input_folder, '*.png'));

for i = 1:length(files)
    img_name = files(i).name;
    img_path = fullfile(input_folder, img_name);
    
    try
        img = imread(img_path);
        
        % convert RGB to grayscale
        img = rgb2gray(img);
        
        [h, w] = size(img);
        if h < target_size(1) || w < target_size(2)
            fprintf('-> [%d/%d]: %s -> Dimension %dx%d \n', i, length(files), img_name, h, w);
            continue;
        end

        fprintf('-> [%d/%d]: %s \n', i, length(files), img_name);
        
        % center crop the image to 1024x1024
        r_start = floor((h - target_size(1)) / 2) + 1;
        c_start = floor((w - target_size(2)) / 2) + 1;
        
        img_cropped = img(r_start:(r_start + target_size(1) - 1), c_start:(c_start + target_size(2) - 1));
        
        [~, base_name, ~] = fileparts(img_name);
        
        % save_path_png = fullfile(output_folder_png, [base_name, '_1024x1024.png']);
        % imwrite(img_cropped, save_path_png);
        
        % save as .mat file to match the original authors format
        I7 = double(img_cropped);
        save_path_mat = fullfile(output_folder_data, [base_name, '.mat']);
        save(save_path_mat, 'I7');
        
    catch ME
        fprintf('Err: %s\n', ME.message);
    end
end
