function frac_out = GS_reflectivity(lambda_in,N_mirror,plotQ)
%
%   Reflecivity of Gemini telescope mirrors vs wavelength.
%   Currently based on data file reflectivity_fresh.dat, from GS M1 2010 ID witness
%   sample, ignoring wavelengths > 1000 nm. Starts at 300 nm.
%
%   Input parameters:
%   -----------------
%   lambda_in   : single value or column vector of wavelength values, in *nm*
%   N_mirror    : Number of mirrors - either 2 or 3 depending on which focus used.
%                                     But allows N_mirror = 1 for testing
%   plotQ       : 1 for plot, 0 for no plot
%
%   Output parameters:
%   ------------------
%   frac_out    : vector, fractional throughput of the mirrors, range 0 - 1, same
%                                                   length as lambda_in
%
%                                                       JGR 17 June 2019
%
%   Presets
%
    lam_min = 300;
    lam_max = 1000;
%
    lambda_in = lambda_in(:).';  % ensures row vector 
    [dim1,~] = size(lambda_in);
    assert(dim1 == 1,'lambda is not a scalar or vector!')
    assert(min(lambda_in)>=lam_min,'lambda value(s) below blue limit!')
    assert(max(lambda_in)<=lam_max,'lambda value(s) above red limit!') 
    assert(N_mirror == 1 || N_mirror == 2 || N_mirror == 3,'Invalid number of mirrors!')
%    
    load('RefData','Reflec_fresh') 
    Reflec_data = Reflec_fresh;      % allow simple change of input file name later ...
%
    Reflec_data(:,2) = Reflec_data(:,2)/100;     % convert from % to fraction
%   
%   Do a spline fit to 1 nm spaced points, over whole range. 
%
    xgrid = (lam_min:lam_max)';  
    spline_data=spline(Reflec_data(:,1),Reflec_data(:,2),xgrid);
%
    if plotQ
        clf
        plot(xgrid,spline_data)
        hold on
        plot(Reflec_data(:,1),Reflec_data(:,2),'kd','MarkerFaceColor','black','MarkerSize', 2, 'DisplayName','input')
    end
%
%   Now use linear interpolation to find values at the desired lambda points
%
    reflec = interp1((lam_min:lam_max),spline_data,lambda_in,'linear');
%
    frac_out = reflec.^N_mirror; 
%
    if plotQ
        #plot(lambda_in, frac_out(:,2), 'color','green', 'LineWidth',2)
        plot(lambda_in, frac_out, 'color', 'green', 'LineWidth', 2, 'DisplayName', sprintf('output with %d mirrors', N_mirror))
        grid on
        xlabel('wavelength (nm)')
        ylabel('reflectivity')
        legend('Location','southeast')
        title('GS\_reflectivity.m')
    end
%    
    return   
end
