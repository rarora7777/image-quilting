function output = synthesis(filename)
    inputRgb = im2double(imread(filename));


    if(length(size(inputRgb)) ~= 3)
        inputRgb = repmat(inputRgb,[1 1 3]);
    end
    inputGray = rgb2gray(inputRgb);
    [m,n] = size(inputGray);

    magnification = 2;
    % size of block is w×w
    w = 50;
    % size of overlap
    o = ceil(w/6);
    % size of output image (w/ additional `o` rows and columns)
    M = ceil(m*magnification/w)*w+o;
    N = ceil(n*magnification/w)*w+o;
    outputGray = zeros(M,N);
    outputRgb = zeros(M,N,3);

    % allInputPatches = extractPatches(inputGray, floor((w+o)/2), inf);


    for i = 1:floor(M/w)
        for j = 1:floor(N/w)
            % mask is 1 for pixels overlapping w/ existing blocks
            mask = zeros(w+o,w+o);
            curOutPatch = outputGray((i-1)*w+1:i*w+o,(j-1)*w+1:j*w+o);

            % Top-left block is picked arbitrarily
            if(i==1 && j==1)
                outputGray(1:w+o,1:w+o) = inputGray(1:w+o,1:w+o);
                outputRgb(1:w+o,1:w+o,:) = inputRgb(1:w+o,1:w+o,:);
                continue;

            elseif(i==1)
                % In first row of blocks, overlaps are only along the left
                % columns
                % x---
                % x---
                % x---
                % x---
                mask(:,1:o) = 1;

                % Get a similar patch by comparing (L2-difference) patches of 
                % the input image with the overlapping portions of the output
                % synthesized so far.
                [nearPatchGray,nearPatchRgb] = getSimilarPatch(...
                    inputRgb, inputGray, ...
                    curOutPatch, mask);

                % For each overlapping pixel, errorVertOverlap contains the 
                % squared-L2 difference between the existing block and the new 
                % block.
                errorVertOverlap = (nearPatchGray.*mask-curOutPatch.*mask).^2;
                errorVertOverlap = errorVertOverlap(:,1:o);

                % `cost` gives the cumulative min. cut cost for each pixel,
                % starting from the bottom. `path` encodes the direction for
                % the min. cost cut.
                [cost,path] = findBoundaryHelper1(errorVertOverlap);
                boundary = zeros(w+o,w+o);
                % Get the position of the topmost pixel of the min. cut
                [~,idx] = min(cost(1,:));
                % And then extract the cut
                boundary(:,1:o) = findBoundaryHelper2(path,idx);

            elseif(j==1)
                % In first column of blocks, overlaps are only along the top
                % rows
                % xxxx
                % ----
                % ----
                % ----
                mask(1:o,:) = 1;

                [nearPatchGray,nearPatchRgb] = getSimilarPatch(...
                    inputRgb, inputGray, ...
                    curOutPatch, mask);

                % Works similar to above.
                errorHorzOverlap = (nearPatchGray.*mask-curOutPatch.*mask).^2;
                errorHorzOverlap = errorHorzOverlap(1:o,:);
                [cost,path] = findBoundaryHelper1(errorHorzOverlap');
                boundary = zeros(w+o,w+o);
                [~,idx] = min(cost(1,:));
                boundary(1:o,:) = (findBoundaryHelper2(path,idx))';

            else
                % In general, top rows as well as left columns overlap
                % xxxx
                % x---
                % x---
                % x---
                mask(:,1:o) = 1;
                mask(1:o,:) = 1;

                [nearPatchGray,nearPatchRgb] = getSimilarPatch(...
                    inputRgb, inputGray, ...
                    curOutPatch, mask);

                errorOverlap = (nearPatchGray.*mask-curOutPatch.*mask).^2;

                errorVertOverlap = errorOverlap(1:o,:);
                [cost1,path1] = findBoundaryHelper1(errorVertOverlap');

                errorHorzOverlap = errorOverlap(:,1:o);
                [cost2,path2] = findBoundaryHelper1(errorHorzOverlap);

                cost = cost1(1:o,:)+cost2(1:o,:);

                boundary = zeros(w+o,w+o);
                [~,idx] = min(diag(cost));
                boundary(1:o,idx:w+o) = (findBoundaryHelper2(path1(idx:o+w,:),o-idx+1))';

                boundary(idx:o+w,1:o) = findBoundaryHelper2(path2(idx:o+w,:),idx);

                boundary(1:idx-1,1:idx-1) = 1;
            end

            % Smoothly merge the two patches at the boundary, if desired.
            % Currently, just choosing the boolean boundary.
            smoothBoundaryGray = imgaussfilt(boundary, 1.5); %boundary;
%             smoothBoundaryGray = .5*ones(size(boundary));
            smoothBoundaryRgb = repmat(smoothBoundaryGray,[1 1 3]);

            newOutPatch = ...
                curOutPatch   .* (smoothBoundaryGray) + ...
                nearPatchGray .* (1-smoothBoundaryGray);

            idxR = (i-1)*w+1:i*w+o;
            idxC = (j-1)*w+1:j*w+o;
            outputGray(idxR, idxC) = newOutPatch;

            outputRgb(idxR, idxC, :) = ...
                outputRgb(idxR, idxC, :) .* (smoothBoundaryRgb) + ...
                nearPatchRgb             .* (1-smoothBoundaryRgb);
        end
    end

    % imshow(a); truesize; figure;
    output = outputRgb(1:M-o,1:N-o,:);
    % imshow(output);truesize;
end
