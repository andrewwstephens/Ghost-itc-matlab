function rate0 = incident_flux(lambda, AB, plotQ)
   log_f_lambda = -AB/2.5 - 2*log10(lambda*10) - 2.408/2.5;  % get f_lambda
   f_lambda = 10.^log_f_lambda;  % in erg /s /cm^2 /Å
   % photon rate incident on top of atmosphere, phot/cm^2/Å/s
   rate0 = f_lambda.*lambda*1E-9/(6.6256E-27*2.99792E8);  
   if plotQ
      plot(lambda, rate0)
      xlabel('wavelength (nm)')
      ylabel('photons / cm^2 / Å / s')
   end
   return   
end
