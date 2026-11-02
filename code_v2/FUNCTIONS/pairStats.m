function st = pairStats(drifter, adcp)
%PAIRSTATS  Paired drifter vs ADCP statistics, NaN-safe.
%   bias and RMSE use (drifter - ADCP):
%   positive bias = drifter faster / larger than ADCP.

drifter = drifter(:);
adcp    = adcp(:);
ok = ~isnan(drifter) & ~isnan(adcp);

st.n = nnz(ok);

if st.n >= 3
    c = corrcoef(drifter(ok), adcp(ok));
    st.r = c(1,2);
else
    st.r = NaN;
end

if st.n > 0
    d = drifter(ok) - adcp(ok);
    st.bias = mean(d);
    st.rmse = sqrt(mean(d.^2));
    st.mean_abs_drifter = mean(abs(drifter(ok)));
    st.mean_abs_adcp    = mean(abs(adcp(ok)));
else
    st.bias = NaN;
    st.rmse = NaN;
    st.mean_abs_drifter = NaN;
    st.mean_abs_adcp    = NaN;
end

end
