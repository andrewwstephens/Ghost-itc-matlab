function RD = nmperpix(lambda,plotQ)
%
%   Function to provide values of wavelength interval (nm) per pixel at any
%   wavelength. Since this varies along each order (due to the anamorphic factor), it
%   is not an overall smooth function of wavelength, and we have to find which order
%   is being used and whereabouts in the order the wanted wavelength is.
%   In this respect it is similar to the blaze function or sum(B^2) calculation, and
%   this routine is based on B2lookup.m. Again, with only a small amount of input
%   data, a 2D interpolation has to be performed to find values at any wavelength.
%   This is done as successive 1D interpolations. The blue and red cameras are done
%   separately.
%                                                                                                     
%   Input parameters:
%   ----------------
%   lambda   : wavelength at which sum(B^2) function is required. Single scalar or 
%                    vector of values, nm.
%
%   Output parameters:
%   ------------------
%   RD       : reciprocal dispersion, nm/pixel. Same size as lambda
%
%                                         G. Robertson  22 July 2019  [ref GHOST 3 141]
%
%
%   Preset
%   
    pix_width = 15;      % pixel width in um
%
%   Reciprocal dispersion results from last table of CST_BGT_007 (for HR, but both
%   same)
% 
    blue_RD =[2.1011E-04	1.8861E-04	1.8247E-04; 1.7097E-04	1.5699E-04	1.5181E-04; ...
        1.4197E-04	1.3265E-04	1.2913E-04];  % these data in nm/um
    blue_order = [64 64 64; 78 78 78; 95 95 95];
    blue_lambda = [533.52 537.72 541.92; 438.38 441.21 444.03; 360.347 362.25 364.16];
%
    red_RD = [4.4521E-04	3.5424E-04	3.4704E-04; 3.1858E-04	2.7323E-04	2.6186E-04; ...
                  2.0791E-04	1.8836E-04	1.8136E-04];
    red_order = [34 34 34; 45 45 45; 66 66 66];
    red_lambda = [997.29 1012.20 1027.10; 756.26 764.76 773.25; 517.47 521.43 525.38];
%    
    [n_triples,~] = size(red_order);  % assumes same data size in red and blue ...
%
%   Get 'phase' of data lambdas, assuming centre of each triple is zero (i.e. order
%   centre). Then get parabolic fit for each row of 3 values (i.e at each order m)
%
    nu_data = zeros(1,3);
    for i = 1:n_triples          
         nu_data = pi*red_order(i,:).*(red_lambda(i,:) - red_lambda(i,2))/red_lambda(i,2);                                            
         poly_vsnu(i,1:3,1) = polyfit(nu_data,red_RD(i,:),2);  % 1st plane for red
    end
    for i = 1:n_triples          
         nu_data = pi*blue_order(i,:).*(blue_lambda(i,:) - blue_lambda(i,2))/blue_lambda(i,2);                                               
         poly_vsnu(i,1:3,2) = polyfit(nu_data,blue_RD(i,:),2);
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
            RD_nu_data(j) = polyval(poly_vsnu(j,:,colour),nu_use(i));
        end
        RD_nu_data = (RD_nu_data(:).')';  % ensures column vector 
%   
%   Now do the interpolation along m axis, using RD evaluated at the correct nu for
%   each data m
%
        switch colour                       % get the interpolated values at required nu. 
            case 1      
                poly_B2_nu = polyfit(red_order(:,2),RD_nu_data,2);   % red
            case 2      
                poly_B2_nu = polyfit(blue_order(:,2),RD_nu_data,2);  % blue
        end
        RD(i) = polyval(poly_B2_nu,m_use(i))*pix_width;
% test
%         m_test = 34:95;
%         RD_test = polyval(poly_B2_nu,m_test);
%         plot(m_test,RD_test)
%         hold on
%         plot(blue_order,RD_nu_data,'kd','MarkerFaceColor','black','MarkerSize',3)
%         plot(red_order,RD_nu_data,'kd','MarkerFaceColor','black','MarkerSize',3)
%         plot(m_use(i),polyval(poly_B2_nu,m_use(i)),'gd','MarkerFaceColor','green','MarkerSize',4)
% end test
    end

    if plotQ
        clf
        plot(lambda, RD, 'LineWidth', 1.5)
        grid on
        xlabel('Wavelength (nm)')
        ylabel('Reciprocal dispersion (nm/pixel)')
        title('nmperpix.m')
        #xlim([363 1000])
        saveas(1, '/tmp/nmperpix.png', 'png');
    end

    return
%
end
