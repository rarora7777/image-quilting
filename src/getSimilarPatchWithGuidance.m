function [nearPatchGray,nearPatchRgb] = getSimilarPatchWithGuidance(...
    alpha,...
    srcRgb, srcGray,...
    tgtGrayPatch, curOutPatch,...
    mask, corrMask, ...
    oldPatch, matchWithOld)

    MAX_ERROR_RELATIVE_THRESHOLD = 1.1;
    [m, n] = size(mask);
    
    % L2-error is given by
    % ∑_i (a_i - b_i)^2 = ∑_i a_i^2 + ∑_i b_i^2 - ∑_i 2*a_i*b_i
    % Compute the three terms in the RHS above:
    outPatchSqrMask = curOutPatch.*curOutPatch.*mask;
    srcGraySqr = srcGray.*srcGray;
    srcGraySqrMask = filter2(mask, srcGraySqr, 'valid');
    outTimesSrcMask = filter2(curOutPatch.*mask, srcGray, 'valid');
    
    overlapMatchError = ...
        (sum(outPatchSqrMask(:)) + srcGraySqrMask - 2*outTimesSrcMask);
    
    
    % Compute the same terms but with srcGrayPatch and tgtGrayPatch and no
    % mask.
    tgtPatchSqrCorr = tgtGrayPatch.*tgtGrayPatch.*corrMask;
    srcPatchSqrCorr = filter2(corrMask, srcGraySqr, 'valid');
    tgtTimesSrcCorr = filter2(tgtGrayPatch.*corrMask, srcGray, 'valid');
    
    
    correspondenceError = ...
        (sum(tgtPatchSqrCorr(:)) + srcPatchSqrCorr - 2*tgtTimesSrcCorr);
    
    if matchWithOld
        oldPatchSqr = oldPatch.*oldPatch;
        oneMat = ones(m, n);
        srcPatchSqrFull = filter2(oneMat, srcGraySqr, 'valid');
        oldTimesSrcFull = filter2(oldPatch, srcGray, 'valid');
        existingMatchError = ...
            sum(oldPatchSqr(:)) + srcPatchSqrFull - 2*oldTimesSrcFull;
        
        overlapMatchError = overlapMatchError + existingMatchError;
    end
    
    
    % err = alpha * overlapMatchError + (1-alpha) * correspondenceError
    errors = ...
        alpha     * overlapMatchError + ...
        (1-alpha) * correspondenceError;
    
    minerror = abs(min(errors(:)));
    
    % select subset of patches where error is below the threshold
    [r, c] = find(errors <= minerror * MAX_ERROR_RELATIVE_THRESHOLD);
    
    % index a random patch
    randint = randi([1 length(r)],1);
    
    
    nearPatchGray = ...
        srcGray(r(randint):r(randint)+m-1, c(randint):c(randint)+n-1);
    nearPatchRgb = ...
        srcRgb(r(randint):r(randint)+m-1, c(randint):c(randint)+n-1,:);
end
