function [sig_xmean,sig_xsig] = SS_calc_sigonmean_sigonsig(lines_mean,lines_std)

nlines = length(lines_mean);

nsim = 1e6;
N    = 26;

for ilines = 1:nlines
    mu0   = lines_mean(ilines);
    sig0  = lines_std(ilines);
    x     = mu0 + sig0*randn(N,nsim);
    xmean = mean(x,1);
    xsig  = std(x,1);
    
    sig_xmean(ilines) = std(xmean);
    sig_xsig(ilines)  = std(xsig); 
end

