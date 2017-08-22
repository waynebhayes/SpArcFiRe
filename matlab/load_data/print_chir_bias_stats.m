function print_chir_bias_stats(chir_given, chir_meas, meas_name)

if any(chir_given == 0)
    error('true chirality measures should not have zero-values\n');
end

% z-values for 95%, 99%, and 99.9% confidence intervals
z950 = norminv(0.975, 0, 1);
z990 = norminv(0.995, 0, 1);
z999 = norminv(0.9995, 0, 1);

chirAgrees = ((chir_given > 0) & (chir_meas > 0)) | ...
    ((chir_given < 0) & (chir_meas < 0));

chirDisagrees = ((chir_given > 0) & (chir_meas < 0)) | ...
    ((chir_given < 0) & (chir_meas > 0));

fprintf('\n');
fprintf('chirality bias-investigation statistics (%s):\n', meas_name);
fprintf('-------------------------------------------------------------\n');

fprintf('distribution of given chiralities:\n');
print_proportions(chir_given, 'given');

fprintf('distribution of measured chiralities:\n');
print_proportions(chir_meas, 'measured');

fprintf('distribution when measure DISAGREES with given (measured values are the opposite):\n');
print_proportions(chir_given(chirDisagrees), 'given');

fprintf('distribution when measure AGREES with given:\n');
print_proportions(chir_given(chirAgrees), 'given');

fprintf('given chirality for instances with indefinite measured chirality:\n');
print_proportions(chir_given(chir_meas == 0), 'given');

fprintf('paired t-test for given vs measured chiralities:\n');
print_paired_ttest(chir_given, chir_meas);

fprintf('paired t-test for disagreeing vs agreeing chiralities:\n');
print_paired_ttest(chir_given(chirAgrees), chir_given(chirDisagrees));

fprintf('-------------------------------------------------------------\n');

function print_proportions(spl_chir, chir_type)
    if isempty(spl_chir)
        fprintf('\tno instances of this type\n');
        return;
    end
    
    n_cw = sum(spl_chir < 0);
    n_ccw = sum(spl_chir > 0);
    n_tot = n_cw + n_ccw;
    pct_cw = (100 * n_cw / n_tot);
    pct_ccw = (100 * n_ccw / n_tot);
    fprintf('\t%d of %d (%2.4f%%) are %s as clockwise\n', ...
        n_cw, n_tot, pct_cw, chir_type);
    fprintf('\t%d of %d (%2.4f%%) are %s as counterclockwise\n', ...
        n_ccw, n_tot, pct_ccw, chir_type);
    p_cw = pct_cw / 100;
    bin_sd = sqrt(n_tot * p_cw * (1 - p_cw));
    fprintf('\tbinomial standard deviation: %2.4f (%2.4f%%)\n', ...
        bin_sd, 100 * bin_sd / n_tot);
    fprintf('\tstandard deviations from 50/50: %2.4f\n', ...
        abs(n_cw - (n_tot/2)) / bin_sd);
    zm = sqrt((p_cw * (1 - p_cw)) / n_tot);
    fprintf('\tconfidence intervals for %% clockwise (using Normal approximation):\n');
    fprintf('\t\tp = 0.05:  [%2.4f%%, %2.4f%%]\n', ...
        pct_cw - 100 * (z950 * zm), pct_cw + 100 * (z950 * zm));
    fprintf('\t\tp = 0.01:  [%2.4f%%, %2.4f%%]\n', ...
        pct_cw - 100 * (z990 * zm), pct_cw + 100 * (z990 * zm));
    fprintf('\t\tp = 0.001: [%2.4f%%, %2.4f%%]\n', ...
        pct_cw - 100 * (z999 * zm), pct_cw + 100 * (z999 * zm));
end

function print_paired_ttest(spl1_chir, spl2_chir)
    spl1_chir = spl1_chir(spl1_chir ~= 0);
    spl2_chir = spl2_chir(spl2_chir ~= 0);
    
    n1 = length(spl1_chir);
    n2 = length(spl2_chir);
    p1 = sum(spl1_chir < 0) / n1;
    p2 = sum(spl2_chir < 0) / n2;
    
    sd = sqrt( ((p1 * (1-p1)) / n1) + ((p2 * (1-p2)) / n2) );
    dp = p1 - p2;
    
    fprintf('\tobserved difference: %2.4f\n', dp);
    fprintf('\tp = 0.05:  [%2.4f, %2.4f]\n', ...
        dp - z950 * sd, dp + z950 * sd);
    fprintf('\tp = 0.01:  [%2.4f, %2.4f]\n', ...
        dp - z990 * sd, dp + z990 * sd);
    fprintf('\tp = 0.001: [%2.4f, %2.4f]\n', ...
        dp - z999 * sd, dp + z999 * sd);
end

end