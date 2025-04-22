function [radius, velocity] = inlet(h_i, h_w, q)
    % Calculate radius and velocity at source
    %
    % Function to calculate inlet velocity and radius from effective pressure
    % Assumes semi-circular cross-section
    % See: Eqn S4, Everett et al. (2018) DOI: 10.1038/s41598-018-31875-8
    %
    % Args:
    %     h_i (float) :   ice thickness (m)
    %     h_w (float) :   water depth (m)
    %     q (float)  :    discharge (m3/s)
    %
    % Returns:
    %     radius (float), velocity (float)
    %
    % Raises:
    %     Error if values produce effective pressure greater than zero
    %

    % Define constants
    RHO_I = 917; % Density of ice (kg/m^3)
    RHO_W = 1000; % Density of water (kg/m^3)
    G = 9.81; % Gravity (m/s^2)
    GLEN_A = 6.e-24; % Flow law creep parameter
    GLEN_N = 3.0; % Glen's flow law exponent
    L = 334.e3; % Latent heat of fusion
    F = 0.1; % Friction factor
    C1 = 1./(RHO_I*L);
    C2 = 2.*GLEN_A*GLEN_N^(-GLEN_N);
    C3 = 2.^0.25*(pi+2.)^0.5/(pi^0.25*(RHO_W*F)^0.5);
    GLEN_N = 3.0; % Glen's flow law exponent

    % Calculate effective pressure
    n_eff = (RHO_I * G * h_i - RHO_W * G * h_w);
    assert(n_eff > 0, 'Terminus is floating - need shallower water or more ice!');
    
    % Calculate the factor
    factor = (C1 / (C2 * C3^2 * n_eff^GLEN_N))^(2 / 7);
    
    % Calculate the area and radius
    area = factor * q^(6 / 7);
    radius = sqrt(2 * area / pi);
    
    % Calculate the velocity
    velocity = q / area;
end