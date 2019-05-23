function r = fSimulate_response_prepare_coef(nGridSpace, nFunctions)
% this function generates nFunctions functions in the spatial domain (nGridSpace)

% pre-define the smoothness
smoothness = 2;

% random function here in spatial domain
r = nan(nFunctions,nGridSpace);
for ii=1:nFunctions
    r(ii,:) = fRandomFunction_Fourier(nGridSpace,smoothness,1);
end

end