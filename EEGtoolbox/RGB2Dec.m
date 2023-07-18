function [H]=RGB2Dec(C)

H=rgb2hex(C);
H=['FF' H(2:end)];
H=hex2dec(H);
