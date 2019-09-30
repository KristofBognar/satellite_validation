function [ instrument_name ] = NDACC_instrument( sitename )
% Provides which type of instrument is taking measurements at the NDACC
% site

switch sitename
    case 'Eureka'
        instrument_name = '125HR';
    case 'NyAlesund'
        instrument_name = '120HR';
    case 'Kiruna'
        instrument_name = '120HR';
    case 'Harestua'
        instrument_name = '120M';
    case 'Pokerflat'
        instrument_name = '120HR';
    case 'Thule'
        instrument_name = '120M';
end

end