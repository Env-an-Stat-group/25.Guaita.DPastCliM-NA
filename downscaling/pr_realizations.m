function Ohat = pr_realizations(EOhat,u_mat,...
    path_main,name_var,suffix,name_model)

Ohat = nan(size(EOhat));

% infer sigma2 from u_mat
sigma2 = var(u_mat,0,2);

% Step 1: compute raw realizations ----------------------------
for i_mth = 1:12
    % load monthly Ot from model file    
    load(fullfile(path_main,['downscaling_models_' name_model],...
        [name_var '_PCR_models_mth=' int2str(i_mth) '-' int2str(i_mth) suffix '.mat']),'O_t_local')
    idx_mth = i_mth:12:size(EOhat,2);
    Ot = O_t_local.O_t;
    Ohat(:,idx_mth) = exp(u_mat(:,idx_mth)-sigma2/2) .* (EOhat(:,idx_mth)+Ot);
end

% Step 2: Hybrid Zero-Inflated Gamma post-processing (per station, per month) ----
[nStation,nTime] = size(Ohat);
Ohat_proc = Ohat;  % start from raw values

for i_mth = 1:12
    idx_mth = i_mth:12:nTime;   % time indices for this calendar month

    for i = 1:nStation
        row = Ohat(i,idx_mth);

        if ~all(isnan(row))
            % Step 2a: negatives → dry
            row(row<0) = 0;

            % Step 2b: define wet-day mask
            pos_vals = row(row>0);
            p = numel(pos_vals)/numel(row);  % wet-day probability

            if numel(pos_vals) >= 10
                % fit Gamma to positives
                phat = gamfit(pos_vals);  % [k,theta]
                k = phat(1); theta = phat(2);
            else
                % fallback: exponential with mean of positives
                if ~isempty(pos_vals)
                    m = mean(pos_vals);
                    k = 1; theta = m;
                else
                    k = 1; theta = 1;  % completely dry station-month
                end
            end

            % Step 2c: resample very high extremes (>99th percentile, and anyway more than 20 mm/day)
            q99 = prctile(pos_vals, 99);
            extreme_mask = row > q99 & row>20;
            row(extreme_mask) = gamrnd(k,theta,[1,sum(extreme_mask)]);

            % Step 2d: optional – keep wet-day values below q99 as-is
        end

        Ohat_proc(i,idx_mth) = row;
    end
end

Ohat = Ohat_proc;

end
