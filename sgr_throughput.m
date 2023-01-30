function frac_out = sgr_throughput(lambda_in,plotQ)
%
%   Throughput of GHOST spectrograph vs wavelength.
%   Based on data from file GHOST_throughput_master.dat from Richard McDermid's PDD
%   ETC files. Using col 1 (wavelength) and col 4 (spectrograph predicted), extracted
%   to Ghost_sgr_tx Matlab array. This can be updated if a better data set becomes
%   available. Wavelengths 360 - 1000 nm. Assumes first and last wavelengths (at
%   least) are integer nm.
%
%   Input parameters:
%   -----------------
%   lambda_in   : single value or column vector of wavelength values, in *nm*
%   plotQ       : 1 for plot, 0 for no plot
%
%   Output parameters:
%   ------------------
%   frac_out    : vector, fractional throughput, range 0 - 1
%
%                                                       JGR 24 June 2019
%
%   Changed to v2 of spectrograph throughput data - 15/1/20
%
%   Get data
%
    load('RefData','Ghost_sgr_tx_v2')
    in_data = Ghost_sgr_tx_v2;      % allow simple change of input file name later ...
    lam_min = min(in_data(:,1));
    lam_max = max(in_data(:,1));
% 
    [dim1,dim2] = size(lambda_in);
    assert(dim1 == 1 | dim2 == 1,'lambda is not a scalar or vector!')
    lambda_in = lambda_in(:).';  % ensures row vector 
    assert(min(lambda_in)>=lam_min,'lambda value(s) below blue limit!')
    assert(max(lambda_in)<=lam_max,'lambda value(s) above red limit!')   
%   
%   Do a spline fit to 1 nm spaced points, over whole range. 
%
    xgrid = (lam_min:lam_max)';  
    spline_data=spline(in_data(:,1),in_data(:,2),xgrid);
%
    if plotQ
        clf
        plot(in_data(:,1),in_data(:,2),'kd','MarkerFaceColor','black','MarkerSize',2)
        hold on
    end
%
%   Now use linear interpolation to find values at the desired lambda points
%
    frac_out = interp1((lam_min:lam_max),spline_data,lambda_in,'linear');
%
    if plotQ
        plot(lambda_in,frac_out, 'color','green', 'LineWidth',2)
        grid on
        xlabel('Wavelength (nm)')
        ylabel('Spectrograph Throughput')
        title('sgr\_throughput.m')
    end
%    
    return   
end
