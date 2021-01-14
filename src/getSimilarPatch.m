function [patchGray, patchRgb] = ...
    getSimilarPatch(inputRgb, inputGray, curOutPatch, mask)
    MAX_ERROR_RELATIVE_THRESHOLD = 1.1;
    
    [m,n] = size(mask);
    
    outPatchSquared = curOutPatch.*curOutPatch.*mask;
    inPatchSquared = inputGray.*inputGray;
    inPatchSquared = filter2(mask, inPatchSquared, 'valid');
    outTimesIn = filter2(curOutPatch.*mask, inputGray, 'valid');
    % L2-error is given by
    % ∑_i (a_i - b_i)^2 = ∑_i a_i^2 + ∑_i b_i^2 - ∑_i 2*a_i*b_i
    errors = sum(outPatchSquared(:)) + inPatchSquared -2*outTimesIn;
    minerror = abs(min(errors(:)));

%     [m,n] = size(mask);
%     m = floor(m/2);
%     n = floor(n/2);
%     
%     curOutPatch = curOutPatch(:)';
%     mask = mask(:)';
%     errors = sum((mask.*(curOutPatch - allInputPatches).^2), 2);
%     minerror = abs(min(errors));
%     errors = reshape(errors, size(inputGray));
    
    % select subset of patches where error is below the threshold
    [x,y] = find(errors <= minerror * MAX_ERROR_RELATIVE_THRESHOLD);
    
    % index a random patch
    randint = randi([1 length(x)],1);
    
    
    % extract the patch
    patchGray = inputGray(x(randint):x(randint)+m-1,y(randint):y(randint)+n-1);
    patchRgb = inputRgb(x(randint):x(randint)+m-1,y(randint):y(randint)+n-1,:);
%     patchGray = inputGray(x(randint)-m:x(randint)+m,y(randint)-n:y(randint)+n);
%     patchRgb = inputRgb(x(randint)-m:x(randint)+m,y(randint)-n:y(randint)+n,:);
end
