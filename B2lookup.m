function B2 = B2lookup(lambda,resoln,flag, plotQ)
%
%   Function to provide values of sum(B^2) at any wavelength - i.e. the (square of)
%   conversion factor from S/N per pixel to S/N per resolution element (for weak
%   spectral features). Since this varies along each order (due to the anamorphic
%   factor), it is not an overall smooth function of wavelength, and we have to find
%   which order is being used and whereabouts in the order the wanted wavelength is.
%       In this respect it is similar to the blaze function calculation, and this
%   routine is based on BlazeFunction.m. But with only a small amount of input
%   data, a 2D interpolation has to be performed to find values at any wavelength.
%   This is done as successive 1D interpolations.
%   The function changes significantly between blue and red cameras, so they have to
%   be done separately.
%   [Uses approximations which assume orders m are all high. Based on Schroeder
%   'Astronomical Optics' p335]         GHOST 3 135.
%                                                                
%   Input parameters:
%   ----------------
%   lambda   : wavelength at which sum(B^2) function is required. Single scalar or 
%                    vector of values, nm.
%   resoln   : 'SR' or 'HR' depending on whether standard or high resolution mode
%   flag     : 1 for interpolated values of sum(B^2), 0 for FWHM
%
%   Output parameters:
%   ------------------
%   B2       : sum of B^2 across the LSF. Same size as lambda
%
%                                                       G. Robertson 19 July 2019
%
%   Added option to produce FWHM instead of sum(B^2). GHOST 3 152, 1 Aug 2019
% 
%
    assert(flag == 0 || flag == 1,'Illegal mode flag value!')
%
%   Results from GHOST 3 18 and FWHM's from the zmx_an2 zemax runs
%
    SR_blue_B2 = [2.5602 2.4959 2.2177; 2.5287 2.4404 2.2449; 2.4524 2.3751 2.2205];
    SR_blue_FWHM = [3.5002 3.3990 3.0481; 3.4505 3.3407 3.0687; 3.3448 3.2432 3.0354];
    SR_blue_order = [64  64  64;  78  78  78;  95  95  95];
    SR_blue_lambda =[541.9300  537.7200  533.5300;  444.0400  441.2100  438.3800; ...
                     364.1600  362.2600  360.350];  %eyeball est
%
    SR_red_B2 = [1.9907 2.4788 2.5355; 2.0782 2.4463 2.5356; 2.1819 2.411 2.5077];
    SR_red_FWHM = [2.7424 3.3808 3.468; 2.8521 3.3427 3.4632; 2.9833 3.2905 3.4204];
    SR_red_order = [36 36 36; 45 45 45; 66 66 66];
    SR_red_lambda = [942.6771 955.9578 969.2385; 756.2718 764.7692 773.2667; ...
                     517.4839 521.4341 525.3843];
%
%   HR results from GHOST 3 12 and 14; FWHMs and lambdas from last table of CST_BGT_007
%
    HR_blue_B2 = [1.6154 1.5650 1.4460; 1.5898 1.5406 1.4197; 1.5399 1.4918 1.4084];
    HR_blue_FWHM = [2.1180 2.0180 1.9180; 2.0580 2.0070 1.8680; 1.9910 1.9330 1.8480];
    HR_blue_order = [64 64 64; 78 78 78; 95 95 95];
    HR_blue_lambda =[541.92 537.72 533.52; 444.03 441.21 438.38; 364.16 362.25 360.347];
%
    HR_red_B2 = [1.3325 1.5631 1.6019; 1.3655 1.5367 1.5979; 1.3915 1.5171 1.5739];
    HR_red_FWHM = [1.772 2.029 2.102; 1.814 2.000 2.091; 1.842 1.968 2.046];
    HR_red_order = [34 34 34; 45 45 45; 66 66 66];
    HR_red_lambda = [997.29 1012.20 1027.10; 756.26 764.76 773.25; 517.47 521.43 525.38];
%
%   Select input data
%
    switch resoln
        case 'SR'
            m_data_red = SR_red_order; 
            m_data_blue = SR_blue_order;
            lambda_data_red = SR_red_lambda; 
            lambda_data_blue = SR_blue_lambda;
            B2_data_red = SR_red_B2*flag + SR_red_FWHM*~flag; 
            B2_data_blue = SR_blue_B2*flag + SR_blue_FWHM*~flag;
        case 'HR'
            m_data_red = HR_red_order; 
            m_data_blue = HR_blue_order;
            lambda_data_red = HR_red_lambda; 
            lambda_data_blue = HR_blue_lambda;
            B2_data_red = HR_red_B2*flag + HR_red_FWHM*~flag; 
            B2_data_blue = HR_blue_B2*flag + HR_blue_FWHM*~flag;
        otherwise
            error('Invalid resolution category input to B^2 interpolation routine!')
    end
%    
    [n_triples,~] = size(m_data_red);  % assumes same data size in red and blue ...
%
%   Get 'phase' of data lambdas, assuming centre of each triple is zero (i.e. order
%   centre). Then get parabolic fit for each row of 3 values (i.e at each order m)
%
    nu_data = zeros(1,3);
    for i = 1:n_triples          
         nu_data = pi*m_data_red(i,:).*(lambda_data_red(i,:) - lambda_data_red(i,2))...
                                                            /lambda_data_red(i,2);
         poly_vsnu(i,1:3,1) = polyfit(nu_data,B2_data_red(i,:),2);  % 1st plane for red
    end
    for i = 1:n_triples          
         nu_data = pi*m_data_blue(i,:).*(lambda_data_blue(i,:) - lambda_data_blue(i,2))...
                                                            /lambda_data_blue(i,2);
         poly_vsnu(i,1:3,2) = polyfit(nu_data,B2_data_blue(i,:),2);
    end
%
%   Get echelle order data for the input lambdas 
%
    load('RefData','Echelle_orders')
    m_vec = Echelle_orders(:,1);
    lambda_B_vec = Echelle_orders(:,2);
    [m_len,~] = size(Echelle_orders);
    lambda = lambda(:).';  % ensures row vector 
    [dim1,n_lambda] = size(lambda);
    assert(dim1 == 1,'lambda is not a scalar or vector!')
%
%   To find which lambda_B is closest to each element of lambda, start by making
%   matrix of the differences. Then use min to find closest.
%
    for i = 1:m_len
        Diff_matrix(i,:) = abs(lambda - lambda_B_vec(i));
    end   
    [~,I] = min(Diff_matrix);  % I is row vector with same length as lambda, giving 
%                                element in Echelle_orders which has closest value
%
    lambda_B_use = lambda_B_vec(I);       % lambda_B_use has same length as lambda
    m_use = m_vec(I);                     % ditto m_use
    lambda_B_use = lambda_B_use(:).';     % Change them back to row vectors
    m_use = m_use(:).';     
    nu_use = pi*m_use.*(lambda - lambda_B_use)./lambda_B_use;  % phase difference wrt order centre                       
%
%   Now do interpolation along order axis at nu of each lambda, to get interpolated value 
%   for each required lambda.
%
    for i = 1:n_lambda
        if m_use(i) >=65, colour = 2; else colour = 1; end  % colour = 1 for red, 2 for blue
        for j = 1:n_triples
            B2_nu_data(j) = polyval(poly_vsnu(j,:,colour),nu_use(i));
        end
        B2_nu_data = (B2_nu_data(:).')';  % ensures column vector 
%   
%   Now do the interpolation along m axis, using B2 evaluated at the correct nu for
%   each data m
%
        switch colour                       % get the interpolated values at required nu. 
            case 1      
                poly_B2_nu = polyfit(m_data_red(:,2),B2_nu_data,2);   % red
            case 2      
                poly_B2_nu = polyfit(m_data_blue(:,2),B2_nu_data,2);  % blue
        end
        B2(i) = polyval(poly_B2_nu,m_use(i));
% test
%         m_test = 34:95;
%         B2_test = polyval(poly_B2_nu,m_test);
%         plot(m_test,B2_test)
%         hold on
%         plot(m_data_blue,B2_nu_data,'kd','MarkerFaceColor','black','MarkerSize',3)
%         plot(m_data_red,B2_nu_data,'kd','MarkerFaceColor','black','MarkerSize',3)
% end test
    end

    if plotQ
       clf
       plot(lambda, B2, 'LineWidth', 1.5)
       grid on
       xlabel('Wavelength (nm)')
       ylabel('S/N per pixel to per res element conversion factor')
       title('B2lookup.m')
       #xlim([363 1000])
       saveas(1, '/tmp/b2lookup.png', 'png');
    end

    return
%
end
