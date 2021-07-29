classdef NanoparticleFunctions
    % functions to simulate nanoparticle polarizability, extinction spectra
    % etc as a function of particle shape, size, geometry
    
    methods
        
        function [hbartau, hbartau1, hbaromega1, hbaromega2, A, B, C, omega_p] = MetalProperties(~, element)
             if strcmp(element,'Au')
                hbartau = 0.071;
                hbartau1 = 0.0716;
                hbaromega1 = 2.43;
                hbaromega2 = 1.52;
                A = 0.132;
                B = -1.755;
                C = 20.43;
                omega_p = 9.06;
            elseif strcmp(element, 'Ag')
                hbartau = 0.021;
                hbartau1 = 0.0760;
                hbaromega1 = 4.02;
                hbaromega2 = 18.5;
                A = -9.71;
                B = -1.111;
                C = 13.77;
                omega_p = 9.17;
            elseif strcmp(element, 'Cu')
                hbartau = 0.103;
                hbartau1 = 0.0528;
                hbaromega1 = 2.12;
                hbaromega2 = 5.43;
                A = -4.36;
                B = -1.655;
                C = 12.31;
                omega_p = 8.88;
            end
        end
        
        function [omegap, f0, Omegap, Gamma0, fj, Gammaj, omegaj] = MetalProperties_LDR(~, element)
             if strcmp(element,'Au')
                omegap = 9.03;
                f0 = 0.760;
                Omegap = (sqrt(f0)*omegap);
                Gamma0 = 0.053;
                fj = ([0.024, 0.010, 0.071, 0.601, 4.384]);
                Gammaj = ([0.241, 0.345, 0.870, 2.494, 2.214]);
                omegaj = ([0.415, 0.830, 2.969, 4.304, 13.32]);
            elseif strcmp(element, 'Ag')
                omegap = 9.01;
                f0 = 0.845;
                Omegap = (sqrt(f0)*omegap);
                Gamma0 = 0.048;
                fj = ([0.065, 0.124, 0.011, 0.840, 5.646]);
                Gammaj = ([3.886, 0.452, 0.065, 0.916, 2.419]);
                omegaj = ([0.816, 4.481, 8.185, 9.083, 20.29]); 
            elseif strcmp(element, 'Cu')
                omegap = 10.83;
                f0 = 0.575;
                Omegap = (sqrt(f0)*omegap);
                Gamma0 = 0.030;
                fj = ([0.061, 0.104, 0.723, 0.638]);
                Gammaj = ([0.378, 1.056, 3.213, 4.305]);
                omegaj = ([0.291, 2.957, 5.300, 11.18]);
             elseif strcmp(element, 'Al')
                omegap = 14.98;
                f0 = 0.523;
                Omegap = (sqrt(f0)*omegap);
                Gamma0 = 0.047;
                fj = ([0.227, 0.050, 0.166, 0.030]);
                Gammaj = ([0.333, 0.312, 1.351, 3.382]);
                omegaj = ([0.162, 1.544, 1.808, 3.473]);               
            end
        end
        
        function [eps_long, V1, a12_long, a14_long, Vparticle, ...
                eps_trans1, eps_trans2, eps_trans3, Vtrans1, ...
                Vtrans2, Vtrans3, a12_trans1, a12_trans2, a12_trans3, ...
                a14_trans1, a14_trans2, a14_trans3] = GeometryFactors(~, shape, R)
             if strcmp(shape, 'Rod')
                % longitudinal mode for Rod
                eps_long = (-1.73.*(R.^(1.45))) - 0.296;
                V1 = 0.896;
                a12_long = 6.92/(1 - eps_long);
                a14_long = (-11/(R.^(2.49))) - 0.0868;

                % transverse modes for Rod
                eps_trans1 = -1.75 + (3.19./(R.^(6.14)));
                eps_trans2 = -0.978 - (0.661./(R.^(1.1)));
                eps_trans3 = -1.57 + (0.0446.*R);
                Vtrans1 = 0.0679 + (1.83./(R.^(2.1)));
                Vtrans2 = 0.891 - (2.28./(R.^(2.53)));
                Vtrans3 = -0.0346 + (0.0111.*R);
                a12_trans1 = 0.0148 + (3.69./(R.^(2.86)));
                a12_trans2 = -21.7 + (22.7./(R.^(0.0232)));
                a12_trans3 = -0.0117 + (0.773/(R.^(1.46)));
                a14_trans1 = 0.0142 - (16.9./(R.^(3.58)));
                a14_trans2 = 1.48 - (3.67./(R.^(0.458)));
                a14_trans3 = -0.256 + (0.0554.*(R.^(0.758)));
                Vparticle = @(R, L) ((pi.*((3.*R) - 1))/(12.*(R.^(3)))).*(L.^3);
            elseif strcmp(shape, 'Triangle')
                % mode for Triangle
                eps_long = (-0.87.*(R.^(1.12))) - 4.33;
                V1 = (-0.645.*(R.^(-1.24))) + 0.678;
                a12_long = 5.57/(1 - eps_long);
                a14_long = (-6.83)/(1 - eps_long);
                Vparticle = @(R, L) ((-0.00544./(R.^2)) + (0.433./R)).*(L.^3);
            elseif strcmp(shape, 'Cage')
                % mode for Cage
                eps_long = (-0.0678.*(R.^(2.02))) - 3.42;
                V1 = (-0.008.*(R.^(2))) + (0.103.*R) + 0.316;
                a12_long = (-0.00405).*(R.^(2.59)) + 2.21;
                a14_long = -13.9;
                Vparticle = @(R, L) (((8.04)./(R.^(3))) - (12./(R.^2)) + (6./R) - 0.00138).*(L.^3);
            elseif strcmp(shape, 'Ellipsoid')
                % mode for Ellipsoid
                eps_long = -0.871 - (1.35.*(R.^(1.54)));
                V1 = 0.994;
                a12_long = 5.52/(1 - eps_long);
                a14_long = (-9.75)./(R.^(2.53));
                Vparticle = @(R, L) (pi./(6.*(R.^2))).*(L.^3);
            elseif strcmp(shape, 'Bicone')
                % mode for Bicone
                eps_long = -0.687 - (2.54.*(R.^(1.5)));
                V1 = (0.648 - 0.441./(R.^(0.687)));
                a12_long = 1.34./(1 - eps_long);
                a14_long = (-1.04)./(1 - eps_long);
                Vparticle = @(R, L) (0.262./(R.^2)).*(L.^3);
            elseif strcmp(shape, 'Disk')
                % mode for Disk
                eps_long = -0.479 - (1.36.*(R.^(0.872)));
                V1 = 0.944;
                a12_long = 7.05./(1 - eps_long);
                a14_long = -10.9./(R.^(0.98));
                Vparticle = @(R, L) ((pi.*(4 + 3.*(R - 1).*(2.*R + pi - 2)))./(24.*(R.^(3)))).*(L.^3);
            elseif strcmp(shape, 'Ring')
               % mode for Ring
               eps_long = 1.39 - (1.31.*(R.^(1.73)));
               V1 = 0.514 + (2.07./(R.^(2.67)));
               a12_long = 7.24./(1 - eps_long);
               a14_long = -19.1./(1 - eps_long);
               Vparticle = @(R, L) (((pi.^2).*(R - 1))./(4.*(R.^(3)))).*(L.^3);
            elseif strcmp(shape, 'Bipyramid')
                % mode for Bipyramid
                eps_long = 1.43 - (4.52.*(R.^(1.12)));
                V1 = 1.96 - (1.73./(R.^(0.207)));
                a12_long = 2.89./(1 - eps_long);
                a14_long = -1.79./(1 - eps_long);
                Vparticle = @(R, L) (0.219./(R.^(2))).*(L.^3);
            elseif strcmp(shape, 'Squared Rod')
                eps_long = -2.28 - (1.47.*(R.^(1.49)));
                V1 = 0.904 - (0.411./(R.^(2.26)));
                a12_long = -0.573 + (3.31./(R.^(0.747)));
                a14_long = 0.213 - (13.1./(R.^(1.97)));
                Vparticle = @(R, L) (1./(R.^(2))).*(L.^3);
            elseif strcmp(shape, 'Cylinder')
                eps_long = -1.59 - (1.96.*(R.^(1.4)));
                V1 = 0.883 - (0.149./(R.^(3.97)));
                a12_long = -1.05 + (3.02./(R.^(0.494)));
                a14_long = 0.0796 - (9.08./(R.^(2.08)));
                Vparticle = @(R, L) (pi./(4.*(R.^(2)))).*(L.^3);
             end
             
             if ~strcmp(shape, 'Rod')
                eps_trans1 = 0;
                eps_trans2 = 0;
                eps_trans3 = 0;
                Vtrans1 = 0;
                Vtrans2 = 0;
                Vtrans3 = 0;
                a12_trans1 = 0;
                a12_trans2 = 0;
                a12_trans3 = 0;
                a14_trans1 = 0;
                a14_trans2 = 0;
                a14_trans3 = 0;
            end

        end
        
        function eps_m = eps_m_gen(~, wavelength, T, cv)
            SF = SolventParameters();
            rhow = SF.WaterDensity(T-273.15);
            n0w = SF.WaterRefractiveIndex(wavelength,T,rhow);
            n0g = SF.GlycerolRefractiveIndex(wavelength);
            n0g(isnan(n0g))=0;
            n_m = sqrt(cv.*n0g.^2 + (1-cv).*n0w.^2);
            eps_m = n_m.^2; % Individual relative permittivies of medium
        end
        
        function eps_p = eps_p_ldr(obj, omega, element)
            % Outputs epsilon of particle according to Lorentz-Drude model parameterised as in
            % Rakić, A. D.; Djurišić, A. B.; Elazar, J. M.; Majewski, M. L.
            % Appl. Opt. 1998, 37 (22), 5271–5283
            % ------ INPUTS ------ %
            % omega (in eV)
            % ------ OUTPUTS ------ %
            % eps_p epsilon of metal at specified wavelengths'
            [omegap, ~, Omegap, Gamma0, fj, Gammaj, omegaj] = obj.MetalProperties_LDR(element);
            eps_r_f = 1 - ((Omegap.^2)./(omega.*(omega.*(1i.*Gamma0))));
            eps_r_b = ((fj.*omegap).^2)./(((omegaj.^2) - (omega'.^2)) + 1i.*omega'.*Gammaj);
            eps_r_b = sum(eps_r_b, 2)';
            eps_p = eps_r_f + eps_r_b;
        end

        
        function eps_p = eps_p_yu(obj, omega, element)
            % Outputs epsilon of particle according to model parameterised as in
            % Yu, R.; Liz-Marzán, L. M.; Abajo, F. J. G. de.
            % Universal Analytical Modeling of Plasmonic Nanoparticles.
            % Chem. Soc. Rev. 2017, 46 (22), 6710–6724.
            % ------ INPUTS ------ %
            % omega (in eV)
            % ------ OUTPUTS ------ %
            % eps_p epsilon of metal at specified wavelengths
            [hbartau, hbartau1, hbaromega1, hbaromega2, A, B, C, omega_p] = obj.MetalProperties(element);            
            eps_b = A + B.*log((hbaromega1 - omega - (1i.*hbartau1))./(hbaromega1 + omega ...
                + (1i.*hbartau1))) + C.*exp(-(omega./hbaromega2));

            eps_p = eps_b - ((omega_p).^2 ./ (omega.*(omega + 1i.*hbartau)));% permittivity of the particle medium
        end
        
        function alpha_j = alphaval(~, s, L, R, a12, a14, V, eps_j, eps_particle, eps_m)
            Aj = @(s, L, R, a12, a14) a12.*(s.^(2)) + ((4*(pi.^2).*(1i.*(V)))./(3.*(L.^(3)))).*(s.^3) ... 
                + a14.*(s.^(4));
            alpha_j = (eps_m./(4.*pi)) .* V .* ...
                (1./((eps_particle./eps_m) - 1) - (1./(eps_j - 1)) - Aj(s, L, R, a12, a14)).^(-1);
        end
        
        function [alpha3, alpha1] = YuPolarizability(obj,wavelength,T,cv,L,R,element,AHFactor,shape)
            % polarizability equations from 
            % Yu, R.; Liz-Marzán, L. M.; Abajo, F. J. G. de.
            % Universal Analytical Modeling of Plasmonic Nanoparticles.
            % Chem. Soc. Rev. 2017, 46 (22), 6710–6724. https://doi.org/10.1039/C6CS00919K.
            % wavelength in nm
            % T is temperature (in Kelvin)
            % cv is volume fraction of glycerol
            % L is length of nanoparticle (in nm)
            % R is L/Width Ratio as defined in Yu et al
            % element is element (Au, Ag, Cu) under consideration
            % AHFactor is ad-hoc factor that sufficienctly replicates
            % experimental prominence of band related to alpha1
            % polarizability tensor element
            % shape is shape to be considered, can be:
            % 'Rod', 'Triangle', 'Cage', 'Ellipsoid', 'Bicone', 'Disk',
            % 'Ring', 'Bipyramid', 'Squared Rod', 'Cylinder'
            c = physconst('lightspeed');
            heV = 4.135667696E-15;
            
            eps_m = obj.eps_m_gen(wavelength, T, cv);
            s = (sqrt(eps_m).*L) ./ wavelength;
            
            energy = ((heV.*c)./(wavelength.*1E-9));
            if strcmp(element, 'Al')
                eps_particle = obj.eps_p_ldr(energy, element);
            else
                eps_particle = obj.eps_p_yu(energy, element);
            end
            
            [eps_long, V1, a12_long, a14_long, Vparticle, ...
                eps_trans1, eps_trans2, eps_trans3, Vtrans1, ...
                Vtrans2, Vtrans3, a12_trans1, a12_trans2, a12_trans3, ...
                a14_trans1, a14_trans2, a14_trans3] = obj.GeometryFactors(shape, R);
            
            alpha3 = obj.alphaval(s, L, R, a12_long, a14_long, Vparticle(R, L).*V1, ...
                eps_long, eps_particle, eps_m);
            
            % correct for only Rod having transverse plasmons
            if strcmp(shape, 'Rod')
                transverseplasmon1 = obj.alphaval(s, L, R, a12_trans1, a14_trans1, Vparticle(R, L).*Vtrans1, ...
                eps_trans1, eps_particle, eps_m);
                transverseplasmon2 = obj.alphaval(s, L, R, a12_trans2, a14_trans2, Vparticle(R, L).*Vtrans2, ...
                eps_trans2, eps_particle, eps_m);
                transverseplasmon3 = obj.alphaval(s, L, R, a12_trans3, a14_trans3, Vparticle(R, L).*Vtrans3, ...
                eps_trans3, eps_particle, eps_m);

                alpha1 = AHFactor*(transverseplasmon1 + transverseplasmon2 + transverseplasmon3);
            else
                alpha1 = zeros(length(wavelength));
            end
            
        end
        
        function [exttotal, extlongitudinal, exttransvers, scattertotal, scatterlongitudinal, scattertransvers...
                , absorptiontotal, absorptionlongitudinal, absorptiontransvers]...
                = YuSpectra(obj,wavelength,shape,T,cv,L,R,element,AHFactor)
            % polarizability equations from 
            % Yu, R.; Liz-Marzán, L. M.; Abajo, F. J. G. de.
            % Universal Analytical Modeling of Plasmonic Nanoparticles.
            % Chem. Soc. Rev. 2017, 46 (22), 6710–6724. https://doi.org/10.1039/C6CS00919K.
            % this function returns the extinction, scattering and
            % absorption spectra
            % wavelength in nm
            % shape is shape to be considered, can be:
            % 'Rod', 'Triangle', 'Cage', 'Ellipsoid', 'Bicone', 'Disk',
            % 'Ring', 'Bipyramid', 'Squared Rod', 'Cylinder'
            % T is temperature (in Kelvin)
            % cv is volume fraction of glycerol
            % L is length of nanoparticle (in nm)
            % R is L/Width Ratio as defined in Yu et al
            % element is element (Au, Ag, Cu) under consideration
            % AHFactor is ad-hoc factor that sufficienctly replicates
            % experimental prominence of band related to alpha1
            % polarizability tensor element            
            eps_m = obj.eps_m_gen(wavelength, T, cv);

            [alpha3, alpha1] = obj.YuPolarizability(wavelength,T,cv,L,R,element,AHFactor,shape);
            
            [~, ~, ~, ~, Vparticle, ...
                ~, ~, ~, ~, ...
                ~, ~, ~, ~, ~, ...
                ~, ~, ~] = obj.GeometryFactors(shape, R);
            
            extinctionprefactor = (8.*(pi.^2)) ./ (sqrt(eps_m) .*wavelength);
            totalplasmons = sum([alpha3; alpha1], 1);
            totaltransverse = 2.*(alpha1);
            exttransvers = (extinctionprefactor .* (imag(totaltransverse)))./Vparticle(R,L);
            extlongitudinal = (extinctionprefactor .* (imag(alpha3)))./Vparticle(R,L);
            exttotal = (extinctionprefactor .* (imag(totalplasmons)))./Vparticle(R,L);
            
            scatteringprefactor = (128 .* (pi.^5))./(3.* (wavelength.^(4)));
            totalplasmons = sum([alpha3; alpha1], 1);
            totaltransverse = 2.*(alpha1);
            scattertransvers = (scatteringprefactor .* (abs(totaltransverse).^2))./Vparticle(R,L);
            scatterlongitudinal = (scatteringprefactor .* (abs(alpha3).^2))./Vparticle(R,L);
            scattertotal = (scatteringprefactor .* (abs(totalplasmons).^2))./Vparticle(R,L);
            
            % absorption cross spectrum calculated from sigma_(abs) =
            % sigma_(ext) - sigma_(sca)
            absorptiontotal = exttotal- scattertotal;
            absorptiontransvers = exttransvers - scattertransvers;
            absorptionlongitudinal = extlongitudinal - scatterlongitudinal;


        end

        function [a1, a3a1, a3] = Yuavalues(obj,wavelength,IntensityWeights,T,cv,L,R,element,AHFactor,shape)
            % polarizability equations from 
            % Yu, R.; Liz-Marzán, L. M.; Abajo, F. J. G. de.
            % Universal Analytical Modeling of Plasmonic Nanoparticles.
            % Chem. Soc. Rev. 2017, 46 (22), 6710–6724. https://doi.org/10.1039/C6CS00919K.
            % wavelength in nm
            % T is temperature (in Kelvin)
            % cv is volume fraction of glycerol
            % L is length of nanoparticle (in nm)
            % R is L/Width Ratio as defined in Yu et al
            % element is element (Au, Ag, Cu) under consideration
            % AHFactor is ad-hoc factor that sufficienctly replicates
            % experimental prominence of band related to alpha1
            % polarizability tensor element
            % shape is shape to be considered, can be:
            % 'Rod', 'Triangle', 'Cage', 'Ellipsoid', 'Bicone', 'Disk',
            % 'Ring', 'Bipyramid', 'Squared Rod', 'Cylinder'
            [alpha3, alpha1] = obj.YuPolarizability(wavelength,T,cv,L,R,element,AHFactor,shape);
            
            % correct for only Rod having transverse plasmons
            if strcmp(shape, 'Rod')
                a3a1 = IntensityWeights(wavelength).*((alpha1.*conj(alpha3)) + ...
                    (conj(alpha1).*alpha3));
                a1 = IntensityWeights(wavelength).*(abs(alpha1).^2);
                a3 = IntensityWeights(wavelength).*(abs(alpha3).^2);

            else
                alpha1 = zeros(length(wavelength));
                a3a1 = IntensityWeights(wavelength).*((alpha1.*conj(alpha3)) + ...
                    (conj(alpha1).*alpha3));
                a1 = IntensityWeights(wavelength).*(abs(alpha1).^2);
                a3 = IntensityWeights(wavelength).*(abs(alpha3).^2);
            end
            
        end
        
        function a3 = Yua3(obj,wavelength,IntensityWeights,T,cv,L,R,element,AHFactor,shape)
            % polarizability equations from 
            % Yu, R.; Liz-Marzán, L. M.; Abajo, F. J. G. de.
            % Universal Analytical Modeling of Plasmonic Nanoparticles.
            % Chem. Soc. Rev. 2017, 46 (22), 6710–6724. https://doi.org/10.1039/C6CS00919K.
            % wavelength in nm
            % T is temperature (in Kelvin)
            % cv is volume fraction of glycerol
            % L is length of nanoparticle (in nm)
            % R is L/Width Ratio as defined in Yu et al
            % element is element (Au, Ag, Cu) under consideration
            % AHFactor is ad-hoc factor that sufficienctly replicates
            % experimental prominence of band related to alpha1
            % polarizability tensor element
            % shape is shape to be considered, can be:
            % 'Rod', 'Triangle', 'Cage', 'Ellipsoid', 'Bicone', 'Disk',
            % 'Ring', 'Bipyramid', 'Squared Rod', 'Cylinder'
            [alpha3, ~] = obj.YuPolarizability(wavelength,T,cv,L,R,element,AHFactor,shape);
            
            a3 = IntensityWeights(wavelength).*(abs(alpha3).^2);
            
        end

        function a3a1 = Yua3a1(obj,wavelength,IntensityWeights,T,cv,L,R,element,AHFactor,shape)
            % polarizability equations from 
            % Yu, R.; Liz-Marzán, L. M.; Abajo, F. J. G. de.
            % Universal Analytical Modeling of Plasmonic Nanoparticles.
            % Chem. Soc. Rev. 2017, 46 (22), 6710–6724. https://doi.org/10.1039/C6CS00919K.
            % wavelength in nm
            % T is temperature (in Kelvin)
            % cv is volume fraction of glycerol
            % L is length of nanoparticle (in nm)
            % R is L/Width Ratio as defined in Yu et al
            % element is element (Au, Ag, Cu) under consideration
            % AHFactor is ad-hoc factor that sufficienctly replicates
            % experimental prominence of band related to alpha1
            % polarizability tensor element
            % shape is shape to be considered, can be:
            % 'Rod', 'Triangle', 'Cage', 'Ellipsoid', 'Bicone', 'Disk',
            % 'Ring', 'Bipyramid', 'Squared Rod', 'Cylinder'
            [alpha3, alpha1] = obj.YuPolarizability(wavelength,T,cv,L,R,element,AHFactor,shape);
            
            % correct for only Rod having transverse plasmons
            if strcmp(shape, 'Rod')
                a3a1 = IntensityWeights(wavelength).*((alpha1.*conj(alpha3)) + ...
                    (conj(alpha1).*alpha3));

            else
                alpha1 = zeros(length(wavelength));
                a3a1 = IntensityWeights(wavelength).*((alpha1.*conj(alpha3)) + ...
                    (conj(alpha1).*alpha3));
            end
            
        end

        function a1 = Yua1(obj,wavelength,IntensityWeights,T,cv,L,R,element,AHFactor,shape)
            % polarizability equations from 
            % Yu, R.; Liz-Marzán, L. M.; Abajo, F. J. G. de.
            % Universal Analytical Modeling of Plasmonic Nanoparticles.
            % Chem. Soc. Rev. 2017, 46 (22), 6710–6724. https://doi.org/10.1039/C6CS00919K.
            % wavelength in nm
            % T is temperature (in Kelvin)
            % cv is volume fraction of glycerol
            % L is length of nanoparticle (in nm)
            % R is L/Width Ratio as defined in Yu et al
            % element is element (Au, Ag, Cu) under consideration
            % AHFactor is ad-hoc factor that sufficienctly replicates
            % experimental prominence of band related to alpha1
            % polarizability tensor element
            % shape is shape to be considered, can be:
            % 'Rod', 'Triangle', 'Cage', 'Ellipsoid', 'Bicone', 'Disk',
            % 'Ring', 'Bipyramid', 'Squared Rod', 'Cylinder'
            [~, alpha1] = obj.YuPolarizability(wavelength,T,cv,L,R,element,AHFactor,shape);
            
            % correct for only Rod having transverse plasmons
            if strcmp(shape, 'Rod')
                a1 = IntensityWeights(wavelength).*(abs(alpha1).^2);

            else
                alpha1 = zeros(length(wavelength));
                a1 = IntensityWeights(wavelength).*(abs(alpha1).^2);
            end
            
        end
       
        function [alpha3, alpha1] = GansPolarizability(obj, element, wavelength, T, cv, xi, AHFactor)
            c = physconst('lightspeed');
            heV = 4.135667696E-15;
            
            
            eps_m = obj.eps_m_gen(wavelength, T, cv);
            energy = ((heV.*c)./(wavelength.*1E-9));
            
            if strcmp(element, 'Al')
                eps_particle = obj.eps_p_ldr(energy, element);
            else
                eps_particle = obj.eps_p_yu(energy, element);
            end


            if xi > 0 % if an ellipsoid i.e. xi > 0
                e = sqrt(1 - ((1+xi).^-2));
                L3_p = ((1-e.^2)./e.^2).*((1./(2.*e)).*log((1+e)./(1-e))-1); % L1 shape parameter for prolate spheroid
                L1_p = (1-L3_p)./2;

                alpha1 = AHFactor.*((eps_particle - eps_m)./(eps_m + L1_p.*(eps_particle - eps_m)));
                alpha3 = (eps_particle - eps_m)./(eps_m + L3_p.*(eps_particle - eps_m));
            else % if a sphere i.e. xi = 0
                alpha3 = (eps_particle - eps_m)./((eps_particle + 2.*(eps_m))./3);
                alpha1 = zeros(1, length(wavelength));
            end

        end
        
        function [exttotal, extlongitudinal, exttransvers, scattertotal, scatterlongitudinal, scattertransvers...
                , absorptiontotal, absorptionlongitudinal, absorptiontransvers]...
                = GansSpectra(obj,wavelength,T,cv,xi,element,AHFactor)
            % wavelength in nm
            % T is temperature (in Kelvin)
            % cv is volume fraction of glycerol
            % xi is shape
            % element is element (Au, Ag, Cu) under consideration
            % AHFactor is ad-hoc factor that sufficienctly replicates
            % experimental prominence of band related to alpha1
            % polarizability tensor element
            eps_m = obj.eps_m_gen(wavelength, T, cv);

            [alpha3, alpha1] = obj.GansPolarizability(element, wavelength, T, cv, xi, AHFactor);
                        
            extinctionprefactor = (8.*(pi.^2)) ./ (sqrt(eps_m) .*wavelength);
            totalplasmons = sum([alpha3; alpha1; alpha1], 1);
            totaltransverse = 2.*(alpha1);
            exttransvers = (extinctionprefactor .* (imag(totaltransverse)));
            extlongitudinal = (extinctionprefactor .* (imag(alpha3)));
            exttotal = (extinctionprefactor .* (imag(totalplasmons)));
            
            scatteringprefactor = (128 .* (pi.^5))./(3.* (wavelength.^(4)));
            totalplasmons = sum([alpha3; alpha1; alpha1], 1);
            totaltransverse = 2.*(alpha1);
            scattertransvers = (scatteringprefactor .* (abs(totaltransverse).^2));
            scatterlongitudinal = (scatteringprefactor .* (abs(alpha3).^2));
            scattertotal = (scatteringprefactor .* (abs(totalplasmons).^2));
            
            % absorption cross spectrum calculated from sigma_(abs) =
            % sigma_(ext) - sigma_(sca)
            absorptiontotal = exttotal- scattertotal;
            absorptiontransvers = exttransvers - scattertransvers;
            absorptionlongitudinal = extlongitudinal - scatterlongitudinal;


        end

        
    end
end

