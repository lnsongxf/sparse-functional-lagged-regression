function rmse = fTransferError(transfer,transfer_true)

% identify lags with non-zero
nonzero_transfers = zeros(2*transfer_true.numOfLags+1,1);
for lag_index = 1:(2*transfer_true.numOfLags+1)
    if max(abs(transfer_true.inSpace(lag_index,:))) > 1e-6 % if it's numerically nonzero, save it as a nonzero transfer functional
        nonzero_transfers(lag_index) = 1;
    end
end

% prepare the structure for RMSE of nonzero functionals
numOfFunctionals = sum(nonzero_transfers);
rmse_transfers = zeros(numOfFunctionals,1);

% cycle over all nonzero functionals
iFunctional = 1;
for ii = 1:(2*transfer_true.numOfLags+1)
    if nonzero_transfers(ii)
        rmse_transfers(iFunctional) = norm(transfer.ONB(ii,:) - transfer_true.ONB(ii,:)) / norm(transfer_true.ONB(ii,:));
        iFunctional = iFunctional + 1;
    end    
end

rmse = mean(rmse_transfers);

end