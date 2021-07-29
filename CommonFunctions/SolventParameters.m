classdef SolventParameters
    %class that calculates water, glycerol densities and refractive indices
        
    methods
        function [rho] = WaterDensity(~, T)
            % Returns density (in g/mL) of pure, air-saturated water for a given temperature (in Celsius) using
            % equation from Jones et al. J. Res. Natl. Inst. Stand. Technol. 97, 335
            % (1992). Valid over the range 5 Celsius <= T <= 40 Celsius.

            rho = 999.84847 + 6.337563e-2*T - 8.523829e-3*T^2 + 6.943248e-5*T^3 - 3.821216e-7*T^4;
            rho = rho/1000;
        end
        
        function [n] = WaterRefractiveIndex(~, lambda, T, rho)
            % This function takes a wavelength (in nanometers), temperature (in
            % Kelvin), and density (in g/mL) and returns the refractive index of
            % water as calculated using the equation in IAPWS "Release on the
            % Refractive Index of Ordinary Water Substance as a Function of Wavelength,
            % Temperature, and Pressure" (1997).
            %
            % The recommended ranges for this function are:
            %
            % -12 Celsius <= T <= 500 Celsius
            % 0 kg*m^-3 <= rho <= 1060 kg*m^-3
            % 200 nm <= lambda <= 1100 nm

            T_star = 273.15;
            rho_star = 1000;
            lambda_star = 0.589;

            T_bar = T/T_star;
            rho_bar = (rho*1000)/rho_star;
            lambda_bar = (lambda/1000)/lambda_star;

            a0 = 0.244257733;
            a1 = 9.74634476e-3;
            a2 = -3.73234996e-3;
            a3 = 2.68678472e-4;
            a4 = 1.58920570e-3;
            a5 = 2.45934259e-3;
            a6 = 0.900704920;
            a7 = -1.66626219e-2;
            lambda_UV = 0.2292020;
            lambda_IR = 5.432937;

            a = 1./rho_bar;
            b = a0 + a1*rho_bar + a2*T_bar + a3*lambda_bar.^2.*T_bar + (a4./lambda_bar.^2) + (a5./(lambda_bar.^2 - lambda_UV^2)) + (a6./(lambda_bar.^2 - lambda_IR^2)) + a7*rho_bar.^2;

            n = sqrt((2*b+a)./(a-b));

        end
        
        function [rho] = GlycerolDensity(~,T)
            % Returns, using linear interpolation, density (in g/mL) of pure glycerol for a given temperature (in Celsius) using
            % data from "Physical Properties of Glycerine and its Solutions" Glycerine Producers' Association: New
            % York, 1963. Valid over the range 0 Celsius <= T <= 290 Celsius.

            % Temperatures
            T_list = [0 10 15 20 30 40 54 75.5 99.5 110 120 130 140 160 180 200 220 240 260 280 290]';
            rho_list = [1.27269 1.26699 1.26443 1.26134 1.25512 1.24896 1.2397 1.2256 1.2097 1.20178...
                1.19446 1.18721 1.17951 1.16440 1.14864 1.13178 1.11493 1.09857 1.08268 1.06725 1.05969]';

            rho = interp1(T_list,rho_list,T);

        end

        function [n] = GlycerolRefractiveIndex(~, lambda)
            % This function takes the data from Birkhoff, R. D. et al. J. Chem. Phys.
            % 69, 4185 (1978) and calulates, using linear interpolation, the refractive
            % index of glycerol for a given wavelength (in nanometers) over the range
            % 51.2 nm - 619.9 nm
            %
            % Above 619.9 nm, the function from Rheims, J. et al. Meas.
            % Sci. Technol. 8, 601 (1997) is used instead, valid over the range >619.9
            % nm - 1050 nm.

            % Light energies in eV as given by Birkhoff et al.
            energies = [2:24 24.2]';
            n_list = [1.47 1.48 1.50 1.56 1.65 1.76 1.87 1.98 2.01 1.94 1.83 1.70 1.57 1.44 1.32 1.22 1.14 1.07 1.01 0.968 0.944 0.939 0.937 0.937];
            k_list = [0 0 0 0 0 0 0 0.121 0.335 0.515 0.646 0.737 0.793 0.809 0.800 0.774 0.736 0.690 0.635 0.576 0.517 0.469 0.437 0.433];

            % Convert to wavelengths in nanometers
            e = 1.6021766208e-19;
            c = 299792458;
            h = 6.62607004e-34;
            wavelengths = (h*c)./(energies*e);
            wavelengths = wavelengths*1e9;
            wavelengths = [1050; wavelengths]; % Pad to range of Rheims data so that conditional function below works properly
            n_list = [1.47 n_list]; % Dummy value to match new range of wavelengths

            n_real_B = @(lambda) interp1(wavelengths,n_list,lambda);
            %k_imag_B = @(lambda) interp1(wavelengths,k_list,lambda);

            A = 1.45797;
            B = 0.00598e-6;
            C = -0.00036e-12;

            % Note that the Rheims curve is adjusted so as to avoid a 0.8%
            % discontinuity with the Birkhoff data at lambda = 619.9 nm
            n_real_R = @(lambda) A + (B./lambda.^2) + (C./lambda.^4) + (n_real_B(619.9) - (A + (B./619.9.^2) + (C./619.9.^4)));

            % fun = @(lambda) (lambda <= 619.9).*(n_real_B(lambda) + k_imag_B(lambda)*1i) + (lambda > 619.9).*n_real_R(lambda);
            fun = @(lambda) (lambda <= 619.9).*n_real_B(lambda) + (lambda > 619.9).*n_real_R(lambda);
            n = fun(lambda);

        end

    end
end

