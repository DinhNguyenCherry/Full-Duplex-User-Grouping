function [varphi, chi, varpi] = GetParameter_Downlink(h, w, phi)


varphi = -log(1 - ( (real(h'*w)^2) ) / ( phi^2 + real(h'*w)^2 ) ) - real(h'*w)^2/(phi^2);

chi = 2*real(h'*w)/(phi^2);

varpi = real( (h'*w)^2 ) / ( phi^2 * ( phi^2 + real(h'*w)^2 ) );

end