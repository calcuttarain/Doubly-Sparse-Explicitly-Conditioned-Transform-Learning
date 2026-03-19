function [settings, datasets_by_n] = generate_settings(T0_list, T1_list, patch_size_list, datasets_input_folders, train_sets, test_sets)
    
    % generate settings
    total_configs = length(datasets_input_folders) * length(T0_list) * length(T1_list) * length(patch_size_list);
    settings(total_configs).T0 = []; 
    
    iter = 0;
    for it0 = 1:length(T0_list)
        for it1 = 1:length(T1_list)
            for ip = 1:length(patch_size_list)
                for ids = 1:length(datasets_input_folders)
                    iter = iter + 1;
                    settings(iter).T0 = T0_list(it0);
                    settings(iter).T1 = T1_list(it1);
                    settings(iter).patch_size = patch_size_list(ip);
                    settings(iter).patch_size_idx = ip;
                    settings(iter).dataset_idx = ids; 
                end
            end
        end
    end
    
    
    % add datasets
    
    num_datasets = length(datasets_input_folders);
    
    datasets_by_n = struct();
    
    for idx_n = 1:length(patch_size_list)
        n = patch_size_list(idx_n);
        
        YH_train_all = cell(1, num_datasets);
        YH_test_all = cell(1, num_datasets);
        
        for d = 1:num_datasets
            folder = datasets_input_folders{d};
            train = train_sets{d};
            test = test_sets{d};
            
            blocks_tr = [];
            blocks_te = [];
            
            for i = 1:length(train)
                filepath = fullfile(folder, [train{i}, '.mat']);
                data = struct2cell(load(filepath));
                img = data{1};
                blocks_tr = [blocks_tr, my_im2col(img, [sqrt(n), sqrt(n)], sqrt(n))];
            end
    
            for i = 1:length(test)
                filepath = fullfile(folder, [test{i}, '.mat']);
                data = struct2cell(load(filepath));
                img = data{1};
                blocks_te = [blocks_te, my_im2col(img, [sqrt(n), sqrt(n)], sqrt(n))];
            end
            
            br_tr = mean(blocks_tr); 
            br_te = mean(blocks_te);
            
            YH_train_all{d} = blocks_tr - (ones(n, 1) * br_tr); 
            YH_test_all{d} = blocks_te - (ones(n, 1) * br_te);
        end
        
        datasets_by_n(idx_n).n = n;
        datasets_by_n(idx_n).YH_train = YH_train_all;
        datasets_by_n(idx_n).YH_test = YH_test_all;
        
    end
end
