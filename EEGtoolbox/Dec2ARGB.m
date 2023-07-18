function [Dec]=Dec2ARGB(Dec)

H=dec2hex(Dec);
if length(H)==8
    H=H(3:end);
end

Dec  = hex2rgb(['#' H],255);