function load_data(path)

    data = load(path);

    if isfield(data, 'paramsin')
        fields = fieldnames(data.paramsin);
        for i = 1:length(fields)
            assignin('base', fields{i}, data.paramsin.(fields{i}));
        end
        
        if isfield(data.paramsin, 'rng_state')
            rng(data.paramsin.rng_state);
        end
    end
    
    if isfield(data, 'paramsout')
        assignin('base', 'paramsout', data.paramsout);
        % fields = fieldnames(data.paramsout);
        % for i = 1:length(fields)
        %     assignin('base', fields{i}, data.paramsout.(fields{i}));
        % end
    end

    if isfield(data, 'ME')
        assignin('base', 'ME', data.ME);
        fprintf('Test Error: %s\n', data.ME.message);
    end

end