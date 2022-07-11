function [varargout] = ghost_etc(varargin)
%
% File generated by IDL2Matlab 1.6 130501 %%

%Initialization of parameters
  I2Mkwn=char('action', 'sn', 'flux', 'texp', 'ab', 'wave', 'nexp', 'seeing', 'airmass', 'mode', 'errfrac', 'nspec', 'nspat', 'quiet', 'help', 'throughput', 'site', 'bright', 'port', 'I2M_pos');
  I2Mkwv={'action', 'sn', 'flux', 'texp', 'ab', 'wave', 'nexp', 'seeing', 'airmass', 'mode', 'errfrac', 'nspec', 'nspat', 'quiet', 'M2I_help', 'throughput', 'site', 'bright', 'port', 'I2M_pos'};
  action=[]; sn=[]; flux=[]; texp=[]; ab=[]; wave=[]; nexp=[]; seeing=[]; airmass=[]; mode=[]; errfrac=[]; nspec=[]; nspat=[]; quiet=[]; M2I_help=[]; throughput=[]; site=[]; bright=[]; port=[]; I2M_pos=[];
  I2M_lst={}; I2M_out=''; lv=length(varargin); if rem(lv,2) ~= 0, I2M_ok=0; else, I2M_ok=1;
  for I2M=1:2:lv; I2M_tmp=varargin{I2M}; if ~ischar(I2M_tmp); I2M_ok=0; break; end; I2Mx=strmatch(I2M_tmp,I2Mkwn); if length(I2Mx) ~=1; I2M_ok=0; break; end; eval([I2Mkwv{I2Mx} '=varargin{I2M+1};']); I2M_lst{(I2M+1)/2}=I2Mkwv{I2Mx}; end; end;
  if ~I2M_ok; for I2M=1:lv; eval([I2Mkwv{I2M} '=varargin{I2M};']); end; end;
  if ~isempty(I2M_pos); for I2M=1:length(I2M_pos); I2Ms=num2str(I2M); I2M_out=[I2M_out 'varargout{' I2Ms '}=' I2M_lst{I2M_pos(I2M)} '; ']; end; end;

%End of parameters initialization

  % Creation of undeclared variables of functions parameters
  fa=1; lext=1; tottpbudg=1; tott=1; snc=1; errstat=1; texpc=1; fluxc=1; 


  %
  % title
  % ghost_etc
  %
  % description
  % routine to compute, in a flexible way, the signal-to-noise ratio (s/n) of a ghost
  % observation of a point continuum source, given different seeing conditions, wavelengths, integration
  % times and so on.
  %
  % inputs
  % action  - determine s/n, exposure time or object flux (s/t/f default s)
  % sn      - input/output s/n as appropriate
  % flux    - input/output object flux in erg/s/cm2 (by a if continuum source) as appropriate
  % texp    - integration time in sec of one exposure (default 1 hour)
  % ab      - if set, then flux is assumed to be in ab magnitudes
  % wave    - reference wavelength (in a) at which to perform calculation (defaut i band 7900 a)
  % nexp    - number of exposure (defaul 1)
  % seeing  - seeing fwhm for a moffat (beta=4) profile
  % airmass - observation airmass (defaulted to 1)
  % mode    - ghost mode (sr, sf, bs, bf, hr, hf, prv)
  % errfrac - output percentage variance contributions from: object, sky, read noise and dark current
  % nspec   - number of spectral pixel to sum up (defaulted to 1), only in continuum
  % nspat   - number of spatial pixel to sum up (defaulted to 1), only in extended
  % quiet   - if set, no print and no plots
  % site    - use extinction data for gemini south ("gs") or north ("gn")
  % bright  - make rough estimate of increased continuum background due to full moon
  % port    - select port on the iss. 1=upward looking, 0=side facing (default). if "s" throughput is selected, the short fiber forces port=1
  %
  % history
  % v1.0 - adapted from muse etc, v1.7 (2007), written by r. mcdermid and
  %        r. bacon.
  % v2.0 - a few important changes:
  %          - corrected the slit length to 75 pixesl in sr/hr modes. 
  %          - changed resolution element to 3.5 pix sr, 0.6x3.5 pix hr
  %          - include option for different sky background
  % v3.0 - removed unused code for other source types
  %      - crude scattered light prescription is now included at top level
  % v4.0 - restructured the code for a cleaner distribution via github.
  %      - added 'very faint' modes to incorporate 1x4 and 1x8 modes
  %      - updated gemini reflectivity to gemini spreadsheet data
  % v5.0 - updated assumed telescope aperture based on engineering data
  %        from gemin (tom hayward)
  %      - updated to use fdr throughput predictsio fom ross, which
  %        include some meaused elements (detector qe, fbpi fibre
  %        measured throughput)
  % v5.1 - added 'port' keyword and 's' throughput to explore shorter fiber and
  %        fixed upward port installation.
  % v5.2 - updated with new throughput data using measured echelle performance.
  %      - removed 'predicted' performance
  %
  M2I_version = '5.2 - 11/10/16';
  %
  %========================================================================================
  %
  % principle common block of variables
  %
  global lbda l lmu extinct tghost fsky flag_ab hc gem ds rn dc daobj dasky M2I_disp  %% share_ghost

  %                          ################################
  %                          # set defaults/keywords/inputs #
  %                          ################################
  % version
  if ( ~keyword_set(quiet))
    [M2I_version] = printt('I2M_a1', 'GHOST_ETC.pro version', 'I2M_a2', M2I_version, 'I2M_pos', [2]);
    % parameter description
  end%if
  if (keyword_set(M2I_help))

    printt('Usage GHOST_ETC, action=S|F|T, flux=flux, sn=sn, texp=texp, /ab, wave=7900, nexp=1,');
    printt('               Seeing=0.8, airmass=1, mode=WFM|NFM, nspec=1, nspat=1');
    if ~isempty(I2M_out),eval(I2M_out);end;return;
  end%if

  if ( ~keyword_set(throughput))
    throughput = 'R';    % assume pdr requirements by default
  end%if
  if ( ~keyword_set(action))
    action = 'S';    % default action = compute s/n
  end%if
  if ( ~keyword_set(flux))
    flux = 1.e-18;    % default flux
  end%if
  if ( ~keyword_set(sn))
    sn = 30.0;    % default s/n
  end%if
  if ( ~keyword_set(texp))
    texp = 3600.0;    % default integration time in sec
  end%if
  if ( ~keyword_set(wave))
    wave = 4500.0;    % default wavelength in angstrom
  end%if
  if ( ~keyword_set(nexp))
    nexp = 1;    % default number of exposure
  end%if
  if ( ~keyword_set(seeing))
    seeing = 0.8;    % default seeing fwhm (moffat, beta=4.0)
  end%if
  if ( ~keyword_set(airmass))
    airmass = 1.0;    % default airmass (zenith)
  end%if
  if ( ~keyword_set(nspec))
    nspec = 1.;    % default number of spectrum pixel to sum up
  end%if
  if ( ~keyword_set(nspat))
    nspat = 1.;    % default number of spatial pixel to sum up
  end%if
  if ( ~keyword_set(mode))
    mode = 'SR';    % default instrument mode
  end%if
  if (keyword_set(ab))
    flag_ab = 1;  
  else
    flag_ab = 0;    % check if ab magnitude flag is set
  end%if
  if ( ~keyword_set(site))
    site = 'GS';    % gs = cerro paranal, gn = mauna kea
  end%if
  if ( ~keyword_set(port))
    port = 0;    % deafult to side port, but can be over-ridden for non-short fiber throughput
  end%if
  if (strcomp(throughput, 'S'))
    port = 1;    % short fiber througput implies the upward port
    % ensure inputs are upper case
  end%if
  mode = strupcase(strtrimi(mode,2));
  site = strupcase(strtrimi(site,2));
  %####################################################################
  % check which mode we are in. this is chosen from:
  %  s - return the s/n for given input source and exposure time
  %  t - return the required exposure time given source and target s/n
  %  f - return the limiting flux for a demand s/n and expo. time
  %####################################################################
  switch action    
    case 'S',

      [texp, flux, wave] = onlyone('I2M_a1', wave, 'I2M_a2', flux, 'I2M_a3', texp, 'I2M_pos', [3, 2, 1]);
      sn = fltarr(eval('n_elements(d1_array(wave,flux,texp))','0') - 2);
      errfrac = fltarr(4,eval('n_elements(sn)','0'));
      if ( ~keyword_set(quiet))
        printt('Compute S/N');
      end%if

    case 'T',

      [sn, flux, wave] = onlyone('I2M_a1', wave, 'I2M_a2', flux, 'I2M_a3', sn, 'I2M_pos', [3, 2, 1]);
      texp = fltarr(eval('n_elements(d1_array(wave,flux,sn))','0') - 2);
      errfrac = fltarr(4,eval('n_elements(texp)','0'));
      if ( ~keyword_set(quiet))
        printt('Compute Exposure time');
      end%if

    case 'F',

      [sn, texp, wave] = onlyone('I2M_a1', wave, 'I2M_a2', texp, 'I2M_a3', sn, 'I2M_pos', [3, 2, 1]);
      flux = fltarr(eval('n_elements(d1_array(wave,texp,sn))','0') - 2);
      errfrac = fltarr(4,eval('n_elements(flux)','0'));
      if ( ~keyword_set(quiet))
        printt('Compute Flux');
      end%if

    otherwise

      [action] = printt('I2M_a1', 'ERROR: Unknown Action (', 'I2M_a2', action, 'I2M_a3', ') should be S, T or F', 'I2M_pos', [2]);
      if ~isempty(I2M_out),eval(I2M_out);end;return;

  end % switch

  %#############################################
  % check if instrument mode parameter is valid
  %#############################################
  switch mode    
    case 'SR',
      if ( ~keyword_set(quiet))
        printt('GHOST Standard Resolution Mode');
      end%if
    case 'SF',
      if ( ~keyword_set(quiet))
        printt('GHOST Standard Resolution Mode - Faint');
      end%if
    case 'SVF',
      if ( ~keyword_set(quiet))
        printt('GHOST Standard Resolution Mode - Very Faint');
      end%if
    case 'BS',
      if ( ~keyword_set(quiet))
        printt('GHOST Beam-Switch Mode');
      end%if
    case 'BF',
      if ( ~keyword_set(quiet))
        printt('GHOST Beam-Switch Mode - Faint');
      end%if
    case 'BVF',
      if ( ~keyword_set(quiet))
        printt('GHOST Beam-Switch Mode - Very Faint');
      end%if
    case 'HR',
      if ( ~keyword_set(quiet))
        printt('GHOST High Resolution Mode');
      end%if
    case 'HF',
      if ( ~keyword_set(quiet))
        printt('GHOST High Resolution Mode - Faint');
      end%if
    case 'PRV',
      if ( ~keyword_set(quiet))
        printt('GHOST PRV mode');
      end%if
    otherwise

      [mode] = printt('I2M_a1', 'Error in parameter Mode: ', 'I2M_a2', mode, 'I2M_pos', [2]);
      if ~isempty(I2M_out),eval(I2M_out);end;return;

  end % switch

  %###################################
  % check if site parameter is valid
  %###################################
  switch site    
    case 'GS',
      if ( ~keyword_set(quiet))
        printt('Assuming Gemini South extinction');
      end%if
    case 'GN',
      if ( ~keyword_set(quiet))
        printt('Assuming Gemini North extinction');
      end%if
    otherwise

      [site] = printt('I2M_a1', 'Error in parameter Site: ', 'I2M_a2', site, 'I2M_pos', [2]);
      if ~isempty(I2M_out),eval(I2M_out);end;return;

  end % switch

  %###################################
  % principle constants and parameters
  %###################################
  % telescope pupil area, in m^2
  %gem = 48.5425  gemini primary useful aperture in m^2
  %gem = 46.0  codr value. no reference
  %gem = 49.1  pdr value. based on optical stops on the primary between 8 and 1.2m
  %gem = !pi*(4.0^2 - 1.0^2)  = 47.1 cdr based on 8m outer diameter and 2m diameter central blockage by baffles+secondary
  gem = 43.748;  % from tom hayward, email to rmcd and
  % mc. this includes the m2 mask, the deployable baffle which has a radius of 1.0 m when
  % in the visible position, and even the obscuration from the spider vanes.
  hc = 6.626075510e-34 .* 299792458.0;  % planck constant * light speed
  pscale = 1.64;  % plate scale arcsec/mm
  %=====================================================================
  % instrument parameters. resolution is as per req.
  %
  % the following parameters are asumed, following analysis by gordon robertson, based on zemax model by john pazder:
  % resolution element = 3.5 pix in sr, 2.2 pix in hr
  % pixels per lens is 4.12 in sr, 2.38 in hr
  standardres = 50000.0;  % standard resolution (sr)
  highres = 75000.0;  % high resolution (hr)
  reselsr = 3.5;  % resolution element in sr (in pixels)
  reselhr = 2.2;  % resolution element in hr (in pixels)
  pixlenssr = 4.12;  % pixels per lens in sr
  pixlenshr = 2.38;  % pixels per lens in hr
  %=====================================================================
  %=====================================================================
  % ghost modes
  % modes are defined in the conops document. main differences are
  % detector binning, spatial sampling and resolution
  switch mode    
    case 'SR',

      res = standardres;      % resolution r=lambda/delta_lambda
      resel = reselsr;      % spectral resolution element in pixels
      nobj = 7.;      % number of object fibers
      nsky = 3.;      % number of sky fibers
      dlens = 240.;      % lens size in micron flat-to-flat
      nslit = nobj .* pixlenssr;      % number of unbinned pixels along object slit
      xbin = 1.;      % detector binning in spectral direction
      ybin = 2.;      % detector binning in spatial direction

    case 'SF',

      res = standardres;      % resolution r=lambda/delta_lambda
      resel = reselsr;      % spectral resolution element in pixels
      nobj = 7.;      % number of object fibers
      nsky = 3.;      % number of sky fibers
      dlens = 240.;      % lens size in micron flat-to-flat
      nslit = nobj .* pixlenssr;      % number of unbinned pixels along object slit
      xbin = 1.;      % detector binning in spectral direction
      ybin = 4.;      % detector binning in spatial direction

    case 'SVF',

      res = standardres;      % resolution r=lambda/delta_lambda
      resel = reselsr;      % spectral resolution element in pixels
      nobj = 7.;      % number of object fibers
      nsky = 3.;      % number of sky fibers
      dlens = 240.;      % lens size in micron flat-to-flat
      nslit = nobj .* pixlenssr;      % number of unbinned pixels along object slit
      xbin = 2.;      % detector binning in spectral direction
      ybin = 4.;      % detector binning in spatial direction

    case 'BS',

      res = standardres;      % resolution r=lambda/delta_lambda
      resel = reselsr;      % spectral resolution element in pixels
      nobj = 7.;      % number of object fibers
      nsky = 7.;      % number of sky fibers
      dlens = 240.;      % lens size in micron flat-to-flat
      nslit = nobj .* pixlenssr;      % number of unbinned pixels along object slit
      xbin = 1.;      % detector binning in spectral direction
      ybin = 2.;      % detector binning in spatial direction

    case 'BF',

      res = standardres;      % resolution r=lambda/delta_lambda
      resel = reselsr;      % spectral resolution element in pixels
      nobj = 7.;      % number of object fibers
      nsky = 7.;      % number of sky fibers
      dlens = 240.;      % lens size in micron flat-to-flat
      nslit = nobj .* pixlenssr;      % number of unbinned pixels along object slit
      xbin = 1.;      % detector binning in spectral direction
      ybin = 8.;      % detector binning in spatial direction

    case 'BVF',

      res = standardres;      % resolution r=lambda/delta_lambda
      resel = reselsr;      % spectral resolution element in pixels
      nobj = 7.;      % number of object fibers
      nsky = 7.;      % number of sky fibers
      dlens = 240.;      % lens size in micron flat-to-flat
      nslit = nobj .* pixlenssr;      % number of unbinned pixels along object slit
      xbin = 2.;      % detector binning in spectral direction
      ybin = 8.;      % detector binning in spatial direction

    case 'HR',

      res = highres;      % resolution r=lambda/delta_lambda
      resel = reselhr;      % spectral resolution element in pixels
      nobj = 19.;      % number of object fibers
      nsky = 7.;      % number of sky fibers
      dlens = 144.;      % lens size in micron, flat-to-flat
      nslit = nobj .* pixlenshr;      % number of unbinned pixels along object slit
      xbin = 1.;      % detector binning in spectral direction
      ybin = 2.;      % detector binning in spatial direction

    case 'HF',

      res = highres;      % resolution r=lambda/delta_lambda
      resel = reselhr;      % spectral resolution element in pixels
      nobj = 19.;      % number of object fibers
      nsky = 7.;      % number of sky fibers
      dlens = 144.;      % lens size in micron, flat-to-flat
      nslit = nobj .* pixlenshr;      % number of unbinned pixels along object slit
      xbin = 1.;      % detector binning in spectral direction
      ybin = 8.;      % detector binning in spatial direction

    case 'PRV',

      res = highres;      % resolution r=lambda/delta_lambda
      resel = reselhr;      % spectral resolution element in pixels
      nobj = 19.;      % number of object fibers
      nsky = 7.;      % number of sky fibers
      dlens = 144.;      % lens size in micron, flat-to-flat
      nslit = nobj .* pixlenshr;      % number of unbinned pixels along object slit
      xbin = 1.;      % detector binning in spectral direction
      ybin = 1.;      % detector binning in spatial direction

  end % switch

  %=====================================================================
  % slit losses
  % compute fraction of flux sampled by the ifu - equivalent to the slit losses
  % this function only considers the geometry of the ifu, no other
  % losses. this should be checked against ross's detailed
  % modelling (but still kept independent from ifu throughput)
  %
  %fracspagauss, mode, seeing, wave, pscale * dlens/1.e3, fa   assumes gaussian seeing
  [fa, wave, seeing, mode] = fracspa('I2M_a1', mode, 'I2M_a2', seeing, 'I2M_a3', wave, 'I2M_a4', pscale .* dlens ./ 1.e3, 'I2M_a5', fa, 'I2M_pos', [5, 3, 2, 1]);  % assumes moffat profile
  %=====================================================================
  % spatial sampling
  %
  dhex = (pscale .* dlens ./ 1.e3) ./ 206265.0;  % diameter of hexagonal lens, flat to flat, in radians
  da = 2. .* sqrt(3.) .* (dhex ./ 2.).^2;  % area of each lens in rad^2
  dasky = nsky .* da;  % total area of sky fibres on sky (rad^2)
  daobj = nobj .* da;  % total area of object fibres on sky (rad^2)
  %=====================================================================
  % spectral sampling
  %
  npix = (resel ./ xbin) .* (nslit ./ ybin);  % total number of pixels summed for object, after binning is applied
  npix_org = resel .* nslit;  % original unbinned pixel number for dark current
  skyfact = 1. + nsky ./ nobj;  % factor for additional pixels used for sky subtraction from simultaneous sky fibers
  nspec = resel;  % number of pixels in resolution element
  ds = (wave .* 1.e-10) ./ resel ./ res;  % size of a spectral pixel in m
  %print,npix,'pixels per object ifu'
  %print,npix*nsky/nobj,'pixels per sky'
  %=====================================================================
  % detector noise properties from requirements
  %
  rn = 4.0;  % readout noise in e-
  dc = 2.6 ./ 3600;  % dark current in e-/sec
  % control amount of output
  switch 1    
    case keyword_set(quiet),
      M2I_disp = 0;
    case eval('n_elements(wave)','0') .* eval('n_elements(flux)','0') .* eval('n_elements(sn)','0') .* eval('n_elements(texp)','0') > 30,

      M2I_disp = 0;
      printt('more than 30 elements computed ... print supressed');

    otherwise
      M2I_disp = 1;
  end % switch

  %                          ###################
  %                          # begin main loop #
  %                          ###################
  n = 0;  % initialize counter for output array
  for i = 0:eval('n_elements(wave)','0') - 1,

    lbda = wave(i +1);    % assign easy variable name
    % check wavelength is ok
    if (((lbda < 3630) | (lbda > 10000)))

      [lbda] = printt('I2M_a1', 'Error Wavelength ', 'I2M_a2', lbda, 'I2M_a3', ' is outside GHOST limits (3630-10000)', 'I2M_pos', [2]);
      if ~isempty(I2M_out),eval(I2M_out);end;return;
    end%if

    l = lbda .* 1.e-10;    % wavelength in m
    lmu = l .* 1.e6;    % wavelength in microns
    %########################################################################
    % atmospheric extinction:
    % interpolate value from the extinction coefficient.
    % gs: patat et al. 2011, a&a, 527, aa91
    % gn: buton, c., copin, y., aldering, g., et al. 2013, a&a, 549, aa8
    %########################################################################
    switch site      
      case 'GS',
        readcol('~/Idl/GHOST_ETC/RefData/paranal_patat11.dat',lext,ext,1);
      case 'GN',
        readcol('~/Idl/GHOST_ETC/RefData/MK_extinction_Buton.dat',lext,ext,'(f,f)',1);
    end % switch

    [extinct, lbda, lext] = interpol('I2M_a1', 10..^(-0.4 .* ext .* airmass), 'I2M_a2', lext, 'I2M_a3', lbda, 'I2M_pos', [3, 2]);
    %########################################################################
    % throughput
    % read throughput data from ross zelhem. this has predicted and budgeted data for cass
    % unit and spectrograph. 10-oct-2016 - includes measured echelle data
    %########################################################################
    readcol('~/Idl/GHOST_ETC/RefData/GHOST_throughput_Echelle.dat',lnm2,tottpbudg,1);
    switch throughput      
      % requirements for cdr
      case 'R',

        if ( ~keyword_set(quiet))
          printt('Assuming throughput specified in Requirement 4110');
        end%if
        switch 1          
          case (lbda >= 3630 & lbda < 3750),
            tghost = 0.08;
          case (lbda >= 3750 & lbda < 4500),
            tghost = 0.127;
          case (lbda == 4500),
            tghost = 0.27;
          case (lbda > 4500 & lbda <= 9000),
            tghost = 0.136;
          case (lbda >= 9000 & lbda <= 9500),
            tghost = 0.055;
          case (lbda > 9500 & lbda <= 10000),
            tghost = 0.021;
          otherwise
            printt('wavelength out of range');
        end % switch


        % total budgeted throughput
      case 'B',

        if ( ~keyword_set(quiet))
          printt('Assuming BUDGETED throughput from Ross Zhelem, 10-10-2016');
        end%if
        [tghost, lbda, tottpbudg] = interpol('I2M_a1', tottpbudg, 'I2M_a2', lnm2 .* 10., 'I2M_a3', lbda, 'I2M_pos', [3, 1]);

      case 'S',

        if ( ~keyword_set(quiet))
          printt('Assuming SHORT FIBER throughput from Mick Edgar, 15-08-2016');
        end%if
        readcol('~/Idl/GHOST_ETC/RefData/GHOSTFullTransmittanceShortFiber.txt',lnm,tott,1);
        [tghost, lbda, tott] = interpol('I2M_a1', tott, 'I2M_a2', lnm .* 10., 'I2M_a3', lbda, 'I2M_pos', [3, 1]);

    end % switch

    %########################################################################
    % gemini mirror reflecitivity
    %
    % various references were used:
    %  codr - flat reflecitivty of 95% for each of three mirrors
    %  pdr - boccas et al. 2006, thin solid films, 502, 275
    %  cdr/fdr - reference data from excel spreadsheet from madeline close, "gemini_reflectivity_2008-2015.txt", version 22 april 2005
    % function for reflectivity of gemini mirrors using data from vucina et al 2008
    % proc. of spie vol. 7012 70122q-1
    % http://www.saao.ac.za/~dod/m5_washing/spie7012_101_gemini.pdf
    %########################################################################
    refdata = 'CDR';    % switch for reflectivity assumption
    %    if refdata eq 1 then m123 = (0.01*(14.09*alog10(lbda) + 36.94))^3 else m123 = 0.925^3
    switch refdata      
      case 'CoDR',
        m123 = 0.95.^3;
      case 'PDR',

        readcol('~/Idl/GHOST_ETC/RefData/boccas_reflectivity.dat',lnmref,reflectivity,1);        % digitized from boccas et al. 2006, thin solid films, 502, 275. used at pdr
        m123 = 0.955 .* (0.01 .* interpol(reflectivity,lnmref .* 10.,lbda)).^3;

      case 'CDR',

        readcol('~/Idl/GHOST_ETC/RefData/GS_M1_Ag_Sample_10-8-2010.dat',lnmref,reflectivity,1,'(f,f)');        % data from excel spreadsheet from madeline close, "gemini_reflectivity_2008-2015.txt", version 22 april 2005
        m1 = 0.95 .* (0.01 .* interpol(reflectivity,lnmref .* 10.,lbda));        % based on worst m1
        m2 = 0.975 .* (0.01 .* interpol(reflectivity,lnmref .* 10.,lbda));        % based on worst m2
        if (port == 1)
          m3 = 1.0;        
        else
          m3 = m1;          % for m3, assume worst m1. skip m3 if upward looking port is set (port=1)
        end%if
        m123 = m1 .* m2 .* m3;        % combined surfaces

    end % switch

    tghost = tghost .* m123;
    %########################################################################
    % sky brightness
    % compute polynomial approximation of no oh paranal sky (in erg/s/a/cm2/arcsec2)
    % the coefficients come from mathcad fit to data from
    % hanuschik r.w., 2003, a&a, 407, 1157
    %########################################################################
    %    bright = 1
    if ( ~keyword_set(quiet))
      if (keyword_set(bright))
        printt('Assuming bright sky');      
      else
        printt('Assuming dark sky conditions');
      end%if
    end%if
    [fsky, bright, lmu] = sky('I2M_a1', lmu, 'bright', bright, 'I2M_pos', [2, 1]);    %,/print)
    fsky = fsky .* 206265.0.^2 .* 1.e7;    % convert in si
    %
    %
    %########################################################################
    % scattered light
    % include a crude scattered light component depending on the
    % wavelength. this scales the effective counts from the source only,
    % and includes this additional flux in the poisson noise calculation
    % (so only adding noise, not signal).
    %########################################################################
    switch 1      
      case lbda <= 3750,
        scatlight = 1.05;
      case lbda > 3750,
        scatlight = 1.02;
    end % switch

    %########################################################################
    % compute requested parameter
    %########################################################################
    switch action      
      case 'S',

        for j = 0:eval('n_elements(flux)','0') - 1,

          for k = 0:eval('n_elements(texp)','0') - 1,

            [scatlight, fa, skyfact, npix_org, npix, nspat, nspec, mode, nexp, errstat, snc] = comp_sn('I2M_a1', snc, 'I2M_a2', flux(j +1), 'I2M_a3', texp(k +1), 'I2M_a4', errstat, 'I2M_a5', nexp, 'I2M_a6', mode, 'I2M_a7', nspec, 'I2M_a8', nspat, 'I2M_a9', npix, 'I2M_b1', npix_org, 'I2M_b2', skyfact, 'I2M_b3', fa, 'I2M_b4', scatlight, 'I2M_pos', [13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 1]);
            sn(n +1) = snc;
            errfrac(:,n +1) = errstat;
            n++;
          end% for

        end% for


      case 'T',

        for j = 0:eval('n_elements(flux)','0') - 1,

          for k = 0:eval('n_elements(sn)','0') - 1,

            [scatlight, fa, skyfact, npix_org, npix, nspat, nspec, mode, nexp, errstat, texpc] = comp_t('I2M_a1', sn(k +1), 'I2M_a2', flux(j +1), 'I2M_a3', texpc, 'I2M_a4', errstat, 'I2M_a5', nexp, 'I2M_a6', mode, 'I2M_a7', nspec, 'I2M_a8', nspat, 'I2M_a9', npix, 'I2M_b1', npix_org, 'I2M_b2', skyfact, 'I2M_b3', fa, 'I2M_b4', scatlight, 'I2M_pos', [13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3]);
            texp(n +1) = texpc;
            errfrac(:,n +1) = errstat;
            n++;
          end% for

        end% for


      case 'F',

        for j = 0:eval('n_elements(texp)','0') - 1,

          for k = 0:eval('n_elements(sn)','0') - 1,

            [scatlight, fa, skyfact, npix_org, npix, nspat, nspec, mode, nexp, errstat, fluxc] = comp_f('I2M_a1', sn(k +1), 'I2M_a2', fluxc, 'I2M_a3', texp(j +1), 'I2M_a4', errstat, 'I2M_a5', nexp, 'I2M_a6', mode, 'I2M_a7', nspec, 'I2M_a8', nspat, 'I2M_a9', npix, 'I2M_b1', npix_org, 'I2M_b2', skyfact, 'I2M_b3', fa, 'I2M_b4', scatlight, 'I2M_pos', [13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 2]);
            flux(n +1) = fluxc;
            errfrac(:,n +1) = errstat;
            n++;
          end% for

        end% for


    end % switch

  end% for

  %                             ##############################
  %                             # plot outputs (if relevant) #
  %                             ##############################
  if (n > 1 &  ~keyword_set(quiet))

    switch action      
      case 'S',

        if ((eval('n_elements(wave)','0') > 1))

          [sn, wave] = plott('I2M_a1', wave, 'I2M_a2', sn, 'xtitle', 'Wavelength', 'ytitle', 'S/N', 'title', 'S/N at Flux: ' + strung(flux) + ' & ' + 'Integ/exp. ' + strung(texp) + 's [ ' + strung(nexp) + ' exp]', 'I2M_pos', [2, 1]);
        end%if

        if ((eval('n_elements(flux)','0') > 1))

          [sn, flux] = plott('I2M_a1', flux, 'I2M_a2', sn, 'xtitle', 'Flux', 'ytitle', 'S/N', 'title', 'S/N at ' + strung(wave) + 'A & Integ/exp. ' + strung(texp) + 's [ ' + strung(nexp) + ' exp]', 'I2M_pos', [2, 1]);
        end%if

        if ((eval('n_elements(texp)','0') > 1))

          [sn, texp] = plott('I2M_a1', texp, 'I2M_a2', sn, 'xtitle', 'Integ/exp. (sec) [' + strung(nexp) + ' exp]', 'ytitle', 'S/N', 'title', 'S/N at ' + strung(wave) + 'A & Flux ' + strung(flux), 'I2M_pos', [2, 1]);
        end%if


      case 'F',

        if ((eval('n_elements(wave)','0') > 1))

          [flux, wave] = plott('I2M_a1', wave, 'I2M_a2', flux, 'xtitle', 'Wavelength', 'ytitle', 'Flux', 'title', 'Flux for S/N ' + strung(sn) + ' & ' + 'Integ/exp. ' + strung(texp) + 's [ ' + strung(nexp) + ' exp]', 'I2M_pos', [2, 1]);
        end%if

        if ((eval('n_elements(sn)','0') > 1))

          [flux, sn] = plott('I2M_a1', sn, 'I2M_a2', flux, 'xtitle', 'S/N', 'ytitle', 'Flux', 'title', 'Flux at ' + strung(wave) + 'A & Integ/exp. ' + strung(texp) + 's [ ' + strung(nexp) + ' exp]', 'I2M_pos', [2, 1]);
        end%if

        if ((eval('n_elements(texp)','0') > 1))

          [flux, texp] = plott('I2M_a1', texp, 'I2M_a2', flux, 'xtitle', 'Integ/exp. (sec) [' + strung(nexp) + ' exp]', 'ytitle', 'Flux', 'title', 'Flux at ' + strung(wave) + 'A & S/N ' + strung(sn), 'I2M_pos', [2, 1]);
        end%if


      case 'T',

        if ((eval('n_elements(wave)','0') > 1))

          [texp, wave] = plott('I2M_a1', wave, 'I2M_a2', texp, 'xtitle', 'Wavelength', 'ytitle', 'Integ/exp. ' + 's [' + strung(nexp) + ']', 'title', 'Int. Time/exp for S/N ' + strung(sn) + ' & Flux ' + strung(flux), 'I2M_pos', [2, 1]);
        end%if

        if ((eval('n_elements(sn)','0') > 1))

          [texp, sn] = plott('I2M_a1', sn, 'I2M_a2', texp, 'xtitle', 'S/N', 'ytitle', 'Integ/exp. ' + 's [' + strung(nexp) + ']', 'title', 'Int. Time/exp at ' + strung(wave) + 'A & Flux ' + strung(flux), 'I2M_pos', [2, 1]);
        end%if

        if ((eval('n_elements(flux)','0') > 1))

          [texp] = plott('I2M_a1', texp, 'I2M_a2', texp, 'xtitle', 'Flux ' + strung(flux), 'ytitle', 'Integ/exp. ' + 's [' + strung(nexp) + ']', 'title', 'Int. Time/exp at ' + strung(wave) + 'A & S/N ' + strung(sn), 'I2M_pos', [2]);
        end%if


    end % switch

  end%if


if ~isempty(I2M_out),eval(I2M_out);end;
 return;
% end of function ghost_etc