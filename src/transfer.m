function output = transfer(textureFile, targetFile)
    srcRgb = im2double(imread(textureFile)); % Input Texture
    tgtRgb = im2double(imread(targetFile)); % Input Image

    if length(size(srcRgb)) ~= 3
        srcRgb = repmat(srcRgb, [1 1 3]);
    end

    if length(size(tgtRgb)) ~= 3
        tgtRgb = repmat(tgtRgb, [1 1 3]);
    end

    srcGray = rgb2gray(srcRgb);
    tgtGray = rgb2gray(tgtRgb);


    [M0, N0] = size(tgtGray);
    

    w = 24;
    o = ceil(w/3);
    alpha = 0.5;    % Only used if NUM_ITERATIONS = 1
    
    NUM_ITERATIONS = 4;

    for iter = 1:NUM_ITERATIONS
        if NUM_ITERATIONS > 1
            alpha = 0.8*(iter-1)/(NUM_ITERATIONS-1) + 0.1;
        end
        
        % M, N, m, n are updated every iteration
        
        % m×n is the number of blocks
        m = ceil((M0-o)/w);
        n = ceil((N0-o)/w);
        
        % M×N is the size of the synthesized image, dependent on the
        % current values of blocksize `w` and overlap width `o`.
        M = m*w + o;
        N = n*w + o;
        
        % Pad `tgtGray` with invalid value (-2) so that the size is M×N.
        temp = -2 * ones(M, N);
        temp(1:M0, 1:N0) = tgtGray(1:M0, 1:N0);
        tgtGray = temp;
        
        % `tgtMask` stores which values are invalid, so that
        % correspondence energy is not computed at those pixels.
        tgtMask = tgtGray > -1;
        
        % In the 1st iteration, simply initialize oldGray, outGray, outRgb
        % to all zeros. Otherwise, retain data from previous iteration but
        % resize to M×N while padding with invalid value (-2).
        if iter > 1
            temp = -2 * ones(M, N);
            prevSize = size(outGray);
            temp(1:prevSize(1), 1:prevSize(2)) = outGray;
            oldGray = zeros(M, N);
            oldGray(1:prevSize(1), 1:prevSize(2)) = outGray;
            outGray = temp;
            temp = -2 * ones(M, N, 3);
            temp(1:prevSize(1), 1:prevSize(2), :) = outRgb;
            outRgb = temp;
        else
            oldGray = zeros(M, N);
            outGray = zeros(M, N);
            outRgb = zeros(M, N, 3);
        end
        
        fprintf("#iter: %d, alpha: %f\n", iter, alpha);
        for i = 1:m
            for j = 1:n
                idxR = (i-1)*w+1:i*w+o;
                idxC = (j-1)*w+1:j*w+o;

                mask = zeros(w+o,w+o);
                curOutPatch = outGray(idxR, idxC);

                if(i==1 && j ==1)
                    [nearPatchGray, nearPatchRgb] = getSimilarPatchWithGuidance(...
                        alpha, ...
                        srcRgb, srcGray, ...
                        tgtGray(idxR, idxC), curOutPatch, ...
                        mask, tgtMask(idxR, idxC), ...
                        oldGray(idxR, idxC), iter>1);
                    outGray(idxR, idxC) = nearPatchGray;
                    outRgb(idxR, idxC,:) = nearPatchRgb;
                    continue;

                elseif(i==1)
                    mask(:,1:o) = 1;
                    [nearPatchGray, nearPatchRgb] = getSimilarPatchWithGuidance(...
                        alpha, ...
                        srcRgb, srcGray,...
                        tgtGray(idxR, idxC), curOutPatch,...
                        mask, tgtMask(idxR, idxC), ...
                        oldGray(idxR, idxC), iter>1);

                    error = (nearPatchGray.*mask-curOutPatch.*mask).^2;
                    error = error(:,1:o);
                    [cost,path] = findBoundaryHelper1(error);
                    boundary = zeros(w+o,w+o);
                    [~,ind] = min(cost(1,:));
                    boundary(:,1:o) = findBoundaryHelper2(path,ind);

                elseif(j==1)
                    mask(1:o,:) = 1;
                    [nearPatchGray, nearPatchRgb] = getSimilarPatchWithGuidance(...
                        alpha, ...
                        srcRgb, srcGray,...
                        tgtGray(idxR, idxC),...
                        curOutPatch,...
                        mask, tgtMask(idxR, idxC), ...
                        oldGray(idxR, idxC), iter>1);

                    error = (nearPatchGray.*mask-curOutPatch.*mask).^2;
                    error = error(1:o,:);
                    [cost,path] = findBoundaryHelper1(error');
                    boundary = zeros(w+o,w+o);
                    [~,ind] = min(cost(1,:));
                    boundary(1:o,:) = (findBoundaryHelper2(path,ind))';

                else
                    mask(:,1:o) = 1;
                    mask(1:o,:) = 1;

                    [nearPatchGray, nearPatchRgb] = getSimilarPatchWithGuidance(...
                        alpha, ...
                        srcRgb, srcGray,...
                        tgtGray(idxR, idxC), curOutPatch,...
                        mask, tgtMask(idxR, idxC), ...
                        oldGray(idxR, idxC), iter>1);

                    error = (nearPatchGray.*mask-curOutPatch.*mask).^2;
                    error1 = error(1:o,:);
                    [cost1,path1] = findBoundaryHelper1(error1');

                    error2 = error(:,1:o);
                    [cost2,path2] = findBoundaryHelper1(error2);

                    cost = cost1(1:o,:)+cost2(1:o,:);

                    boundary = zeros(w+o,w+o);
                    [~,ind] = min(diag(cost));
                    boundary(1:o,ind:w+o) = (findBoundaryHelper2(path1(ind:o+w,:),o-ind+1))';

                    boundary(ind:o+w,1:o) = findBoundaryHelper2(path2(ind:o+w,:),ind);

                    boundary(1:ind-1,1:ind-1) = 1;

                end

%                 smoothBoundaryGray = imgaussfilt(boundary, 1.5, 'Padding', 'replicate');
                smoothBoundaryGray = boundary;
                smoothBoundaryRgb = repmat(smoothBoundaryGray, [1 1 3]);
                
                newOutPatch = ...
                    curOutPatch   .* (smoothBoundaryGray) + ...
                    nearPatchGray .* (1-smoothBoundaryGray);
                
                
                
                outGray(idxR, idxC) = newOutPatch;
                outRgb(idxR, idxC, :) = ...
                    outRgb(idxR, idxC, :) .* (smoothBoundaryRgb) + ...
                    nearPatchRgb          .* (1-smoothBoundaryRgb);
            end
        end

        output = outRgb(1:M0, 1:N0, :);

        w = round(2*w/3);
        o = ceil(w/3);
        
        figure; imshow(output); truesize;
    end
end
