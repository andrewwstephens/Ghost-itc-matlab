function fraction_tx = IFU_trans(Res,FWHM)
%
%   Calculate the fraction transmitted by either the SR or HR IFU, as a function of
%   seeing FWHM. For convenience and speed it uses polynomial fits to curves derived
%   using Monte Carlo integration. The original calculations were in Aperture_thru.m
%   and plotted over a range of seeing using call_Aperture_thru.m
%   The seeing model is a Moffat profile, following formulas from R McDermid in the
%   GHOST PDD. The calculation used the actual multiple-hex shapes of the IFUs.
%
%   Included microlens array throughput factors. SR 0.93, HR 0.85.  9 Jan 2020  
%
%   Input parameters:
%   -----------------
%   Res     : 'SR' for standard resolution, 'HR' for high resolution
%   FWHM    : seeing FWHM in arcsec; single scalar value only
%
%   Output parameters:
%   ------------------
%   fraction_tx : single scalar value of transmitted fraction
%
%   Presets
%   
    poly_SR = [-0.601093757388945,   2.491686724776579,  -3.321623316076719, ...
        0.954991358859481, 0.926565900249881];
    poly_SR_bad = [0.017680338582827,  -0.192379590733923,   0.817695279572324, ...
        -1.662347054268447,  1.461657092371379];
    poly_HR = [0.285444811536887,  -1.867077232744627,   4.552284297394247, ...
        -4.811391795562376, 1.405374986540268,   0.883408649616416];
    poly_HR_bad = [0.011858109451155,  -0.138739933878572,   0.635235598240202, ...
        -1.389283926814856, 1.306482356152507];
    SR_IFU_tx = 0.93;  % SR IFU lenslet throughput factor (Ross, GVR-5124.4)
    HR_IFU_tx = 0.85;  % Equivalent for HR
    
    assert(FWHM>=0,'Seeing FWHM must be non-negative!')
%
    switch Res
        case 'SR'
            if FWHM < 0.28
                 alpha = 0.701*FWHM + 5*eps;  % for best seeing use integration in equiv circle
                fraction_tx = (1 - (1 + (0.333/alpha)^2)^-3);
            elseif FWHM <= 1.55
                fraction_tx = polyval(poly_SR,FWHM);
            elseif FWHM < 3
                fraction_tx = polyval(poly_SR_bad,FWHM);
            else
                alpha = 0.701*FWHM;  % for really bad seeing use integration in equiv circle
                fraction_tx = (1 - (1 + (0.333/alpha)^2)^-3);
            end
            fraction_tx = fraction_tx*SR_IFU_tx;
        case 'HR'
            if FWHM < 0.28
                 alpha = 0.701*FWHM + 5*eps;  
                fraction_tx = (1 - (1 + (0.333/alpha)^2)^-3);
            elseif FWHM <= 1.55
                fraction_tx = polyval(poly_HR,FWHM);
            elseif FWHM < 3
                fraction_tx = polyval(poly_HR_bad,FWHM);
            else
                alpha = 0.701*FWHM;  
                fraction_tx = (1 - (1 + (0.333/alpha)^2)^-3);
            end
            fraction_tx = fraction_tx*HR_IFU_tx;
        otherwise
            error('Invalid resolution type')
    end
    
    fprintf('Fraction of source in the IFU = %0.3f\n', fraction_tx)
    return
%        
end

