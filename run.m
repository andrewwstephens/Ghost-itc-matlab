#!/usr/bin/env octave
SNmain(
   lambda=363:1:1000, % = numpy.arange(363, 1000, 1)
   AB=18,          % AB magnitude of target
   t_exp=900,      % Total exposure time (seconds)
   ZD=30,          % Zenith distance (degrees)
   resoln='SR',    % 'SR' or 'HR'
   seeing=1.0,     % Seeing (arcseconds)
   SB='SB80',      % 'SB20', 'SB50', 'SB80' or 'SBAny'
   N_mirror=3,     % 2 for bottom port, 3 for side port
   SR_sky=3,       % number of SR sky microlenses: 3, 7 or 10
   bin_spat=1,     % spatial binning
   N_sub=1,        % Number of sub-exposures that make up t_exp
   plotQ=1);
