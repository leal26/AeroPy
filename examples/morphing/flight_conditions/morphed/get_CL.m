function [CL]=get_CL(airfoil,alpha,V,H)
% This function will output the coefficient of lift (CL) corresponding to
% the optimal OML at a given set of conditions (i.e., parameter values).
%
% Inputs:
%       airfoil: a string that describes the parent airfoil being
%           considered. Current options are 'NACA0012', 'NACA4415',
%           'NACA641212', and 'glider'. 
%       alpha: a floating point number describing the angle of attack.
%           Valid range is 0-12 degrees. 
%       V: a floating point number describing the airspeed. Valid range is
%           20-65 m/s. 
%       H: a floating point number describing the airspeed. Valid range is
%           5,000-40,000 ft. Default value is 10,000 ft.

if nargin < 4
    H = 10000;
    if nargin < 3
        error('Not enough inputs.')
    end
end

data_calculated = true; % DO NOT CHANGE THIS UNLESS YOU NEED TO RECALCULATE THE DATA (IT TAKES APPROX. 12 HOURS)
if ~data_calculated % if the data is not already calculated:
    % call the subfunction below to calculate the data
    calculate_data; % no inputs/outputs necessary
    % call the subfunction below to generate a kriging function for CL
    generate_kriging_function; % no inputs/outputs necessary
end

filename = append('krig_', airfoil, '.mat');
load(filename, 'model');

% finally calculate
CL = predictor([alpha,V,H], model);

end
