function frac_out = Extinc_Paranal(lambda_in,ZD,plotQ)
%
%   Calculate attenuation due to atmospheric extinction, for Gemini South.
%   Use data from Paranal:  Patat et al AA 527 AA91 2011.
%   The output is fractional throughput, *not* magnitudes.
%
%   Input parameters:
%   -----------------
%   lambda_in   : single value or vector of wavelength values, in *nm*
%   ZD          : zenith distance of the observation, in *degrees*
%   plotQ       : 1 for plot, 0 for no plot
%
%   Output parameters:
%   ------------------
%   frac_out    : row vector giving fractional throughput of the atmosphere, 
%                   range 0 - 1, same length as lambda_in
%
%                                                       JGR 17 June 2019
%   Presets
%
    lam_min = 333;
    lam_max = 1000;
% 
    lambda_in = lambda_in(:).';  % ensures row vector 
    [dim1,~] = size(lambda_in);
    assert(dim1 == 1,'lambda is not a scalar or vector!')
    assert(min(lambda_in)>=lam_min,'lambda value(s) below blue limit!')
    assert(max(lambda_in)<=lam_max,'lambda value(s) above red limit!')
    assert(ZD >= 0 && ZD < 90, 'ZD out of range!')
%
%   Extinction values from Patat et al. These are magnitudes per air mass
%
    lambda = [3325,3375,3425,3475,3525,3575,3625,3675,3725,3775,3825,3875,3925,3975,...
              4025,4075,4125,4175,4225,4275,4325,4375,4425,4475,4525,4575,4625,4675,...
              4725,4775,4825,4875,4925,4975,5025,5075,5125,5175,5225,5275,5325,5375,...
              5425,5475,5525,5575,5625,5675,5725,5775,5825,5875,5925,5975,6025,6075,...
              6125,6175,6225,6275,6325,6375,6425,6475,6525,6575,6625,6675,6725,6775,...
              7060,7450,7940,8500,8675,8850,10000]';
     lambda = lambda/10;  % Convert wavelengths to nm
     k_lam = [0.686,0.606,0.581,0.552,0.526,0.504,0.478,0.456,0.43,0.409,0.386,0.378,...
              0.363,0.345,0.33,0.316,0.298,0.285,0.274,0.265,0.253,0.241,0.229,0.221,...
              0.212,0.204,0.198,0.19,0.185,0.182,0.176,0.169,0.162,0.157,0.156,0.153,...
              0.146,0.143,0.141,0.139,0.139,0.134,0.133,0.131,0.129,0.127,0.128,0.13,...
              0.134,0.132,0.124,0.122,0.125,0.122,0.117,0.115,0.108,0.104,0.102,0.099,...
              0.095,0.092,0.085,0.086,0.083,0.081,0.076,0.072,0.068,0.064,0.064,0.048,...
              0.042,0.032,0.03,0.029,0.022]';       
%
%   Do a spline fit to 1 nm spaced points, over whole range. Do it vs
%   lambda^-4 because more linear so more accurate fit.
%
    xgrid = (lam_min:lam_max)';  
    lambda_mf = lambda.^-4;
    xgrid_mf = xgrid.^-4;
    spline_data=spline(lambda_mf,k_lam,xgrid_mf);
%
    if plotQ
        clf
        plot(xgrid,spline_data)
        hold on
        plot(lambda,k_lam,'kd','MarkerFaceColor','black','MarkerSize',3, 'DisplayName','model')
    end
%
%   Now use linear interpolation to find values at the desired lambda points
%
    extinc = interp1((lam_min:lam_max),spline_data,lambda_in,'linear');
%
    airmass = 1/cosd(ZD);
    frac_out = 10.^(-0.4*airmass*extinc); 
%
    if plotQ
        plot(lambda_in, extinc, 'color', 'red', 'LineWidth', 2, 'MarkerSize', 30, 'DisplayName','interp')
        plot(lambda_in, frac_out, 'color', 'green', 'LineWidth', 2, 'DisplayName', 'Throughput')
        grid on
        xlabel('Wavelength (nm)')
        ylabel('Extinction')
        legend('Location','northwest')
        title('Extinc\_Paranal.m')
        end

    return
end
