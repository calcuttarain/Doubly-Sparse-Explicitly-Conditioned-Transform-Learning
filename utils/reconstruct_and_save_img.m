function reconstruct_and_save_img(W, X, br, patchSize, imageSize, numPatchesPerImage, imageNames, methodName)
    % reconstructie
    W_pinv = pinv(W);
    Y_hat = W_pinv * X;
    
    % adun media
    Y_hat = Y_hat + (ones(size(Y_hat,1),1) * br);
    
    startIdx = 1;
    for i = 1:length(imageNames)
        endIdx = startIdx + numPatchesPerImage - 1;
        Y_hat_img = Y_hat(:, startIdx:endIdx);
        
        I_hat = col2im(Y_hat_img, patchSize, imageSize, 'distinct');
        
        imwrite(uint8(I_hat), sprintf('results/img/%s_%s.png', methodName, imageNames{i}));
        
        startIdx = endIdx + 1;
    end
end