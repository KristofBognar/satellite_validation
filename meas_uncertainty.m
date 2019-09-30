function [ x, e1, e2 ] = meas_uncertainty( m1, m2 )
%MEAS_UNCERTAINTY estimate measurand variability and measurement
%uncertainties from two datasets (Fioletov et al., 2006)
%
% INPUT:
%   coincident measurements of the same quantity from two different instruments
%
% OUTPUT:
%   x: measurand variability
%   e1,e2: measurement uncertainties for both datasets
%
% Kristof Bognar, June 2018


x=sqrt(0.5*abs( var(m1) + var(m2) - var(m1-m2) ));
e1=sqrt(0.5*abs( var(m1) - var(m2) + var(m1-m2) ));
e2=sqrt(0.5*abs( var(m2) - var(m1) + var(m1-m2) ));

end

