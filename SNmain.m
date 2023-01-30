function [SN, vars]=SNmain(lambda,AB,t_exp,ZD,resoln,seeing,SB,N_mirror,SR_sky,bin_spat,N_sub,plotQ)
%
%   Main routine for signal/noise calculation for GHOST.
%   * No provision for spectral binning in this routine *
%
%   Input parameters:
%   -----------------
%   lambda      : wavelength(s) at which S/N is required. Vector or scalar. In *nm*
%   AB          : AB magnitude of object, assumed unresolved
%   t_exp       : exposure time, seconds. *Total* over N_sub exposures if N_sub > 1
%   ZD          : Zenith distance of object, degrees
%   resoln      : Resolution mode - 'SR' for standard, 'HR' for high resolution
%   seeing      : Seeing disc FWHM, arcseconds
%   SB          : Sky brightness class, 'SB20', 'SB50', 'SB80' or 'SBAny'
%   N_mirror    : Number of Gemini reflections: 2 for axial port, 3 for side port
%   SR_sky      : For SR only - number of sky microlenses: 3, 7 or 10
%   bin_spat    : Binning factor in spatial direction; 1 for no binning
%   N_sub       : Number of sub-exposures comprising the total t_exp. 1 for single exp
%   plotQ       : 1 for plot, 0 for no plot
%
%   Output parameters:
%   ------------------
%   SN          : Signal/noise per resolution element at each value of input lambda
%   vars        : array of variances, one value for each lambda. Values in phot^2
%                  row 1: object, row 2 sky, row 3 dark, row 4 readout
%
%                                   G. Robertson  23 July 2019. [GHOST 3 140]
%
%   Presets
%
    #area_tel = 491000.0;  % Gemini primary mirror area, cm^2
    area_tel = 476893.76;   # math.pi * 395**2 - math.pi * 65**2
    
    #RON_red = 2.2;  % red CCD read noise e- rms (NRC rep Jun 2019; slow readout; av of 4)
    #RON_blue = 2.25; % blue CCD read noise e- rms (NRC rep Jun 2019; slow readout; av of 4)
    RON_red = 4.5;  % red CCD read noise e- rms (to match Gemini ITC "Bright Target")
    RON_blue = 4.5; % blue CCD read noise e- rms (to match Gemini ITC "Bright Target")
    
    dark_red = 0.825; % red dark noise, e- /pix /hr (NRC rep Jun 2019; av of 4)
    dark_blue = 1.175; % blue dark noise, e- /pix /hr (NRC rep Jun 2019; av of 4)   
    area_SR = 0.939; % SR IFU area, arcsec^2
    area_HR = 0.917; % HR IFU area, arcsec^2
    lambda_cross = 533.50; % wavelength of crossover from blue to red camera, nm
    Npix_SR = 18.9;  % Number of (unbinned) pixels in length of 1 star SR slit
    Npix_HR = 30.8;  % Number of (unbinned) pixels in length of HR  star slit
%    
    lambda = lambda(:).';  % ensures row vector 
    [dim1,n_lambda] = size(lambda);
    assert(dim1 == 1,'lambda is not a scalar or vector!')   
    assert(min(lambda)>=360,'lambda value(s) below blue limit!')
    assert(max(lambda)<=1000,'lambda value(s) above red limit!')
    assert(N_mirror == 2 || N_mirror == 3,'Illegal number of telescope reflections!')
    assert(abs(round(N_sub) - N_sub)<10*eps && N_sub>0,'Number of subexposures must be integer!')
    vars = zeros(4,n_lambda);
    switch resoln
        case 'SR'
            area_IFU = area_SR;
            slit_len = Npix_SR;
            assert(SR_sky == 3 || SR_sky == 7 || SR_sky == 10,'Illegal number of SR sky lenses!')
            sky_coeff = 1 + 7/SR_sky;
        case 'HR'
            area_IFU = area_HR;
            slit_len = Npix_HR;
            sky_coeff = 3.714;  % i.e 1 + 19/7
        otherwise
            error('Undefined resolution type!')
    end
%
%   Get number of signal photons from the object in each 'pixel' (at lambda values, which
%   won't in general be actual CCD pixel wavelengths, but will use pixel width near that
%   value, which is OK)
%   
    log_f_lambda = -AB/2.5 - 2*log10(lambda*10) - 2.408/2.5;  % get f_lambda
    f_lambda = 10.^log_f_lambda;  % in erg /s /cm^2 /Å
%
    rate0 = f_lambda.*lambda*1E-9/(6.6256E-27*2.99792E8);
               % photon rate incident on top of atmosphere, phot /cm^2 /Å /s
    rate1 = rate0*area_tel;  %     phot /Å /s  for whole telescope area
    extinc_frac = Extinc_Paranal(lambda,ZD,0);  % Paranal extinction; fraction passed
    rate2 = rate1.*extinc_frac; % after atmos extinction; phot /Å /s
    rate3 = rate2.*GS_reflectivity(lambda,N_mirror,0); % after telescope mirrors; phot /Å /s  
    rate4 = rate3.*IFU_trans(resoln,seeing); % drop by IFU injection loss
    rate5 = rate4.*GHOST_cable_eta(lambda);  % drop by Cass Unit and cable throughputs 
    rate6 = rate5.*sgr_throughput(lambda,0); % drop by sgr and detector throughputs    
    rate7 = rate6.*BlazeFunction(lambda,0);    % drop by individual order blaze functions
    RD = nmperpix(lambda,0)*10;  % Reciprocal dispersion, Å/pix
    rate8 = rate7.*RD*0.99;   % rate in phot /s /pixel; drop by 1% loss on slitview beamsplitter
    object = rate8*t_exp;    % phot /pix over whole exposure time; summed over slit length
%
%   Assemble noise variances
%
    sky_flux = Sky_contin(lambda,SB,0);  % sky *continuum* flux, erg /s /Å /cm^2 /arcs^2
    sky_rate1 = sky_flux.*lambda*1E-9*area_tel*t_exp/(6.6256E-27*2.99792E8); % phot /Å /arcs^2
    sky_rate2 = sky_rate1.*GS_reflectivity(lambda,N_mirror,0);   % losses through the system
    sky_rate2a = sky_rate2.*extinc_frac; % Hanuschik sky brightness is for outside atmosphere(!)
    sky_rate3 = sky_rate2a.*GHOST_cable_eta(lambda).*sgr_throughput(lambda,0);
    sky_rate4 = sky_rate3.*BlazeFunction(lambda,0);  
    sky_rate5 = sky_rate4*area_IFU; % phot /Å after injection into IFU
    sky_rate6 = sky_rate5.*RD; % phot /pix over whole exposure time; summed over slit length
    sky_rate7 = sky_rate6*sky_coeff*0.99; % scale up to allow for sky subtraction noise; 1% BS loss
%    
    dark0 = (lambda < lambda_cross)*dark_blue + (lambda >= lambda_cross)*dark_red; %  e- /pix /hr  
    dark1 = dark0*t_exp/3600; %  e- /pix 
    dark2 = dark1*slit_len;   %  e- /pix, summed over slit length
%
    ron_rms = (lambda < lambda_cross)*RON_blue + (lambda >= lambda_cross)*RON_red;
    ron_var = slit_len*N_sub*(ron_rms.^2)/bin_spat; % read noise variance for 1 spectral pixel,
                                                    % all spatial pixels, all subexposures

# Q: Should there be a gain value in here to convert from object photons to electrons?

    var_total = object + sky_rate7 + dark2 + ron_var; % for 1 spectral pixel
    % vars = [object; sky_rate7; dark2; ron_var];

    SN_pix = object./sqrt(var_total);                % S/N per pixel
    SN = SN_pix.*sqrt(B2lookup(lambda,resoln,1,0));  % S/N per resolution element
%       
    if plotQ
        clf
        
        subplot(2,1,1)
        plot(lambda,object, 'color','red', 'LineWidth',1.5, 'DisplayName', 'Source')
        hold on
        plot(lambda,sky_rate7, 'color', 'blue', 'LineWidth',1.5, 'DisplayName', 'Sky')
        plot(lambda,dark2, 'color', 'black', 'LineWidth',1.5, 'DisplayName', 'Dark')
        plot(lambda,ron_var, 'color', 'green', 'LineWidth',1.5, 'DisplayName', 'RN^2')
        grid on
        legend('Location','northwest')
        xlabel('Wavelength (nm)')
        ylabel('electrons (gain=1)')
        xlim([363 1000])
        title('SNmain.m')

        subplot(2,1,2)
        plot(lambda, SN, 'LineWidth',1.5, 'DisplayName', 'S/N / res element')
        hold on
        plot(lambda, SN_pix, 'LineWidth',1.5, 'DisplayName', 'S/N / pixel')
        grid on
        legend('Location','northwest')
        xlabel('Wavelength (nm)')
        ylabel('Signal / Noise')
        xlim([363 1000])
        title(sprintf('AB=%.1f  texp=%.0f  seeing=%.1f  BG=%s  ZD=%d', AB, t_exp, seeing, SB, ZD))
        saveas(1, '/tmp/SN.png', 'png');
    end
return
end
