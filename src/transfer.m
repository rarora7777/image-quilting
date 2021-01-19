function output = transfer(textureFile, targetFile)
    srcRgb = im2double(imread(textureFile)); % Input Texture
    tgtRgb = im2double(imread(targetFile)); % Input Image

    % figure; imshow(a); truesize;
    % figure; imshow(b); truesize;

    if length(size(srcRgb)) ~= 3
        srcRgb = repmat(srcRgb, [1 1 3]);
    end

    if length(size(tgtRgb)) ~= 3
        tgtRgb = repmat(tgtRgb, [1 1 3]);
    end

    srcGray = rgb2gray(srcRgb);
    tgtGray = rgb2gray(tgtRgb);

%     maskSrc = srcGray < -1;
%     maskTgt = tgtGray < -1;
% 
%     srcGray(maskSrc) = -1;
%     tgtGray(maskTgt) = -1;

    [M0, N0] = size(tgtGray);
    
%     outGray = zeros(M0, N0);
%     outRgb = zeros(M0, N0, 3);

    w = 24;
    alpha = 0.5;
    o = round(w/6);
    
    NUM_ITERATION = 4;

    for iter = 1:NUM_ITERATION
        if NUM_ITERATION > 1
            alpha = 0.8*(iter-1)/(NUM_ITERATION-1) + 0.1;
        end
        
        m = ceil((M0-o)/w);
        n = ceil((N0-o)/w);
        M = m*w + o;
        N = n*w + o;
        temp = -2 * ones(M, N);
        temp(1:M0, 1:N0) = tgtGray(1:M0, 1:N0);
        tgtGray = temp;
        tgtMask = tgtGray > -1;
        
        oldGray = zeros(M, N);
        
        if iter > 1
            temp = -2 * ones(M, N);
            prevSize = size(outGray);
            temp(1:prevSize(1), 1:prevSize(2)) = outGray;
            oldGray(1:prevSize(1), 1:prevSize(2)) = outGray;
            outGray = temp;
            temp = -2 * ones(M, N, 3);
            temp(1:prevSize(1), 1:prevSize(2), :) = outRgb;
            outRgb = temp;
        else
            outGray = zeros(M, N);
            outRgb = zeros(M, N, 3);
        end
        
        fprintf("#iter: %d, alpha: %f\n", iter, alpha);
        for i = 1:m
            for j = 1:n
%                 if all(all(maskTgt((i-1)*w+1:i*w+o,(j-1)*w+1:j*w+o)))
%                     outGray((i-1)*w+1:i*w+o,(j-1)*w+1:j*w+o)=0;
%                     outRgb((i-1)*w+1:i*w+o,(j-1)*w+1:j*w+o,:)=0;
%                     continue;
%                 end
                mask = zeros(w+o,w+o);
                curOutPatch = outGray((i-1)*w+1:i*w+o,(j-1)*w+1:j*w+o);

                if(i==1 && j ==1)
                    [nearPatchGray, nearPatchRgb] = getSimilarPatchWithGuidance(...
                        alpha, ...
                        srcRgb, srcGray, ...
                        tgtGray(1:w+o,1:w+o), curOutPatch, ...
                        mask, tgtMask(1:w+o,1:w+o), ...
                        oldGray(1:w+o,1:w+o), iter>1);
                    outGray(1:w+o,1:w+o) = nearPatchGray;
                    outRgb(1:w+o,1:w+o,:) = nearPatchRgb;
                    continue;

                elseif(i==1)
                    mask(:,1:o) = 1;
                    [nearPatchGray, nearPatchRgb] = getSimilarPatchWithGuidance(...
                        alpha, ...
                        srcRgb, srcGray,...
                        tgtGray((i-1)*w+1:i*w+o,(j-1)*w+1:j*w+o), curOutPatch,...
                        mask, tgtMask((i-1)*w+1:i*w+o,(j-1)*w+1:j*w+o), ...
                        oldGray((i-1)*w+1:i*w+o,(j-1)*w+1:j*w+o), iter>1);

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
                        tgtGray((i-1)*w+1:i*w+o,(j-1)*w+1:j*w+o),...
                        curOutPatch,...
                        mask, tgtMask((i-1)*w+1:i*w+o,(j-1)*w+1:j*w+o), ...
                        oldGray((i-1)*w+1:i*w+o,(j-1)*w+1:j*w+o), iter>1);

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
                        tgtGray((i-1)*w+1:i*w+o,(j-1)*w+1:j*w+o), curOutPatch,...
                        mask, tgtMask((i-1)*w+1:i*w+o,(j-1)*w+1:j*w+o), ...
                        oldGray((i-1)*w+1:i*w+o,(j-1)*w+1:j*w+o), iter>1);

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

                smoothBoundaryGray = imgaussfilt(boundary, 1.5);
                smoothBoundaryRgb = repmat(boundary, [1 1 3]);
                
                newOutPatch = ...
                    curOutPatch   .* (smoothBoundaryGray) + ...
                    nearPatchGray .* (1-smoothBoundaryGray);
                
                idxR = (i-1)*w+1:i*w+o;
                idxC = (j-1)*w+1:j*w+o;
                
                outGray(idxR, idxC) = newOutPatch;
                outRgb(idxR, idxC, :) = ...
                    outRgb(idxR, idxC, :) .* (smoothBoundaryRgb) + ...
                    nearPatchRgb          .* (1-smoothBoundaryRgb);
            end
        end

        output = outRgb(1:M0, 1:N0, :);
%         output(repmat(tgtMask, [1 1 3])) = 0;


%         w = round(w*0.7);
        w = round(2*w/3);
        o = round(w/6);
        
%         if NUM_ITERATION > 1
%             alpha = 0.8*(iter-1)/(NUM_ITERATION-1) + 0.1;
%         else
%             continue;
%         end

%         tgtGray = outGray;
%         oldGray = tgtGray;
%         srcGray(maskSrc) = -1;
%         tgtGray(tgtMask) = -1;
%         M = ceil((m-o)/w)*w + o;
%         N = ceil((n-o)/w)*w + o;
%         outGray = zeros(M, N);
%         outRgb = zeros(M, N, 3);
% 
%         temp = -ones(M, N);
%         [mTemp, nTemp] = size(tgtGray);
%         temp(1:mTemp, 1:nTemp) = tgtGray;
%         tgtGray = temp;
%         oldGray = tgtGray;
%         
        figure; imshow(output); truesize;
    end
    
%     output = outRgb(1:M0, 1:N0, :);
end
