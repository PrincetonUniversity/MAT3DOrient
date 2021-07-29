classdef CramerRaoFunctions
    %Cramer Rao lower bound functions for fluorescence and scattering
    %orientation of single particle problem
        
    methods
        function cdist = cosdistance(~, theta1, phi1)
            % cosine error
            cdist = acos(((sin(theta1) .* sin(0)) .* cos(phi1 - 0)) + (cos(theta1) .* cos(0)));
        end

        
        function [n0] = n_m(~,lambda,T,cv)
            % refractive index, T in Kelvin
            SP = SolventParameters;
            rhow = SP.WaterDensity(T-273.15);
            n0w = SP.WaterRefractiveIndex(lambda,T,rhow);
            n0g = SP.GlycerolRefractiveIndex(lambda);
            n0g(isnan(n0g))=0;
            n0 = sqrt(cv.*n0g.^2 + (1-cv).*n0w.^2);
        end
        
        function [Theta, Delta, A, B, C, H] = InstrResp(~,NACond,NAObj,n0)
            % instrument factors from Fourkas & Yang
            Theta = asin(NACond./n0);
            Delta = asin(NAObj./n0);
            A = (1/6)-(1/4)*cos(Delta)+(1/12)*cos(Delta).^3;
            B = (1/8)*cos(Delta) - (1/8)*cos(Delta).^3;
            C = (7/48) - (1/16)*cos(Delta) - ((1/16)*cos(Delta).^2) - ((1/48)*cos(Delta).^3);
            H = cos(2*Theta);
        end
        
        function Jval = Jthetathetaf(~, theta, phi, time, Beta, I, A, B, C, n_detectors)
            % time is time taken to observe photons
            % Beta is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent), should be n_detector x wavelength
            % matrix
            % I is Intensity/(1-SBR^{-1}) (can be wavelength dependent),
            % should be n_detector x wavelength matrix
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            
            if size(I, 2) == 1 % if intensity is a simple scalar
                if length(A) ~= 1 % if we have multiple wavelengths
                    I = I./length(A); % spread intensity across multiple wavelengths
                end
            end
            
            if n_detectors == 3
                I1 = I(1, :); Beta1 = Beta(1, :);
                I2 = I(2, :); Beta2 = Beta(2, :);
                I3 = I(3, :); Beta3 = Beta(3, :);
                I4 = 0; Beta4 = 0;
            else
                I1 = I(1, :); Beta1 = Beta(1, :);
                I2 = I(2, :); Beta2 = Beta(2, :);
                I3 = I(3, :); Beta3 = Beta(3, :);
                I4 = I(4, :); Beta4 = Beta(4, :);
            end
            Term1 = ((A.^2).*(C.^2).*((-1 + Beta1).^2).*I1.*(cos(theta).^2).*(cos(2.*phi).^2).*(sin(theta).^2))...
                ./(Beta1.*((A + B.*(sin(theta).^2)).^3).*(A.*(3 + Beta1) + (B.*(3 + Beta1) ...
                - C.*(Beta1 - 1).*cos(2.*phi)).*(sin(theta).^2)));
            Term2 = ((A.^2).*(C.^2).*((-1 + Beta2).^2).*I2.*(cos(theta).^2).*(cos(2.*phi).^2).*(sin(theta).^2))...
                ./(Beta2.*((A + B.*(sin(theta).^2)).^3).*(A.*(3 + Beta2) + (B.*(3 + Beta2) ...
                + C.*(Beta2 - 1).*cos(2.*phi)).*(sin(theta).^2)));
            Term3 = ((A.^2).*(C.^2).*((-1 + Beta3).^2).*I3.*(sin(2.*phi).^2).*(sin(2.*theta).^2))...
                ./(4.*Beta3.*((A + B.*(sin(theta).^2)).^3).*(A.*(3 + Beta3) + (sin(theta).^2).*(B.*(3 + Beta3)...
                - C.*(Beta3 - 1).*sin(2.*phi))));
            Term4 = ((A.^2).*(C.^2).*((-1 + Beta4).^2).*I4.*(sin(2.*phi).^2).*(sin(2.*theta).^2))...
                ./(4.*Beta4.*((A + B.*(sin(theta).^2)).^3).*(A.*(3 + Beta4) + (sin(theta).^2).*(B.*(3 + Beta4)...
                + C.*(Beta4 - 1).*sin(2.*phi))));
            if n_detectors == 3
                Jval = sum(time.*(Term1 + Term2 + Term3));
            elseif n_detectors == 4
                Jval = sum(time.*(Term1 + Term2 + Term3 + Term4));
            elseif n_detectors == 2
                Jval = sum(time.*(Term1 + Term3));
            end
            Jval = sum(Jval, 2);
        end
        
        function Jval = JthetathetafNB(~, theta, phi, time, I, A, B, C, n_detectors)
            % time is time taken to observe photons
            % here is ideal case - no background
            % I is Intensity (can be wavelength dependent),
            % should be n_detector x wavelength matrix
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            
            if size(I, 2) == 1 % if intensity is a simple scalar
                if length(A) ~= 1 % if we have multiple wavelengths
                    I = I./length(A); % spread intensity across multiple wavelengths
                end
            end
            
            if n_detectors == 3
                I1 = I(1, :);
                I2 = I(2, :);
                I3 = I(3, :);
                I4 = 0;
            else
                I1 = I(1, :);
                I2 = I(2, :);
                I3 = I(3, :);
                I4 = I(4, :);
            end
            
            if n_detectors == 3
                Jval = time.*((1./(4.*(A + B.*(sin(theta).^2)).^3)).*(A.^2).*(C.^2).*...
                    (4.*(cos(theta).^2).*(cos(2.*phi).^2).*((I1./(B + C.*cos(2.*phi) + A.*(csc(theta).^2)))...
                    + (I3./(B - C.*cos(2.*phi) + A.*(csc(theta).^2)))) + ...
                    (sin(2.*theta).^2).*(sin(2.*phi).^2).*((I2./(A + ...
                    (sin(theta).^2).*(B + C.*sin(2.*phi)))))));
            elseif n_detectors == 4
                Jval = time.*((1./(4.*(A + B.*(sin(theta).^2)).^3)).*(A.^2).*(C.^2).*...
                    (4.*(cos(theta).^2).*(cos(2.*phi).^2).*((I1./(B + C.*cos(2.*phi) + A.*(csc(theta).^2)))...
                    + (I3./(B - C.*cos(2.*phi) + A.*(csc(theta).^2)))) + ...
                    (sin(2.*theta).^2).*(sin(2.*phi).^2).*((I2./(A + (sin(theta).^2).*(B + C.*sin(2.*phi))))...
                    + (I4./(A + (sin(theta).^2).*(B - C.*sin(2.*phi)))))));
            end
            Jval = sum(Jval, 2);
        end
        
        
        function Jval = Jthetaphif(~, theta, phi, time, Beta, I, A, B, C, n_detectors)
            % time is time taken to observe photons
            % Beta is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent), should be n_detector x wavelength
            % matrix
            % I is Intensity/(1-SBR^{-1}) (can be wavelength dependent),
            % should be n_detector x wavelength matrix
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            
            if size(I, 2) == 1 % if intensity is a simple scalar
                if length(A) ~= 1 % if we have multiple wavelengths
                    I = I./length(A); % spread intensity across multiple wavelengths
                end
            end
            
            if n_detectors == 3
                I1 = I(1, :); Beta1 = Beta(1, :);
                I2 = I(2, :); Beta2 = Beta(2, :);
                I3 = I(3, :); Beta3 = Beta(3, :);
                I4 = 0; Beta4 = 0;
            else
                I1 = I(1, :); Beta1 = Beta(1, :);
                I2 = I(2, :); Beta2 = Beta(2, :);
                I3 = I(3, :); Beta3 = Beta(3, :);
                I4 = I(4, :); Beta4 = Beta(4, :);
            end
            
            Term1 = -(A.*(C.^2).*((1 - Beta1).^2).*I1.*cot(theta).*(csc(theta).^2).*sin(4.*phi))...
                ./(2.*Beta1.*((B + A.*(csc(theta).^2)).^2).*(-C.*(Beta1 - 1).*cos(2.*phi) + ...
                (3 + Beta1).*(B + A.*(csc(theta).^2))));
            Term2 = -(A.*(C.^2).*((1 - Beta2).^2).*I2.*cot(theta).*(csc(theta).^2).*sin(4.*phi))...
                ./(2.*Beta2.*((B + A.*(csc(theta).^2)).^2).*(C.*(Beta2 - 1).*cos(2.*phi) + ...
                (3 + Beta2).*(B + A.*(csc(theta).^2))));
            Term3 = (A.*(C.^2).*((1 - Beta3).^2).*I3.*cot(theta).*(csc(theta).^2).*sin(4.*phi))...
                ./(2.*Beta3.*((B + A.*(csc(theta).^2)).^2).*((3 + Beta3).*(B + A.*(csc(theta).^2)) ...
                - C.*(Beta3 - 1).*sin(2.*phi)));
            Term4 = (A.*(C.^2).*((1 - Beta4).^2).*I4.*cot(theta).*(csc(theta).^2).*sin(4.*phi))...
                ./(2.*Beta4.*((B + A.*(csc(theta).^2)).^2).*((3 + Beta4).*(B + A.*(csc(theta).^2)) ...
                + C.*(Beta4 - 1).*sin(2.*phi)));
            if n_detectors == 3
                Jval = sum(time.*(Term1 + Term2 + Term3));
            elseif n_detectors == 4
                Jval = sum(time.*(Term1 + Term2 + Term3 + Term4));
            end
            Jval = sum(Jval, 2);
        end
        
        function Jval = JthetaphifNB(~, theta, phi, time, I, A, B, C, n_detectors)
            % time is time taken to observe photons
            % here is ideal case - no background
            % I is Intensity (can be wavelength dependent),
            % should be n_detector x wavelength matrix
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            
            if size(I, 2) == 1 % if intensity is a simple scalar
                if length(A) ~= 1 % if we have multiple wavelengths
                    I = I./length(A); % spread intensity across multiple wavelengths
                end
            end
            
            if n_detectors == 3
                I1 = I(1, :);
                I2 = I(2, :);
                I3 = I(3, :);
                I4 = 0;
            else
                I1 = I(1, :);
                I2 = I(2, :);
                I3 = I(3, :);
                I4 = I(4, :);
            end
            
            if n_detectors == 3
                Jval = time.*((1./(2.*((B + A.*(csc(theta).^2)).^2))).*A.*(C.^2).*cot(theta).*...
                    (csc(theta).^2).*sin(4.*phi).*((-I1./(B + C.*cos(2.*phi) + A.*(csc(theta).^2))) + ...
                    (I2./(B + A.*(csc(theta).^2) + C.*sin(2.*phi))) - ...
                    (I3./(B + A.*(csc(theta).^2) - C.*cos(2.*phi)))));
            elseif n_detectors == 4
                Jval = time.*((1./(2.*((B + A.*(csc(theta).^2)).^2))).*A.*(C.^2).*cot(theta).*...
                    (csc(theta).^2).*sin(4.*phi).*((-I1./(B + C.*cos(2.*phi) + A.*(csc(theta).^2))) + ...
                    (I2./(B + A.*(csc(theta).^2) + C.*sin(2.*phi))) - ...
                    (I3./(B + A.*(csc(theta).^2) - C.*cos(2.*phi))) + ...
                    (I4./(B + A.*(csc(theta).^2) - C.*sin(2.*phi)))));
            end
            Jval = sum(Jval, 2);
        end

        
        function Jval = Jphiphif(~, theta, phi, time, Beta, I, A, B, C, n_detectors)
            % time is time taken to observe photons
            % Beta is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent), should be n_detector x wavelength
            % matrix
            % I is Intensity/(1-SBR^{-1}) (can be wavelength dependent),
            % should be n_detector x wavelength matrix
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            
            if size(I, 2) == 1 % if intensity is a simple scalar
                if length(A) ~= 1 % if we have multiple wavelengths
                    I = I./length(A); % spread intensity across multiple wavelengths
                end
            end
            
            if n_detectors == 3
                I1 = I(1, :); Beta1 = Beta(1, :);
                I2 = I(2, :); Beta2 = Beta(2, :);
                I3 = I(3, :); Beta3 = Beta(3, :);
                I4 = 0; Beta4 = 0;
            else
                I1 = I(1, :); Beta1 = Beta(1, :);
                I2 = I(2, :); Beta2 = Beta(2, :);
                I3 = I(3, :); Beta3 = Beta(3, :);
                I4 = I(4, :); Beta4 = Beta(4, :);
            end

            
            Term1 = ((C.^2).*((Beta1 - 1).^2).*I1.*(sin(2.*phi).^2))...
                ./(Beta1.*(B + A.*(csc(theta).^2)).*(-C.*(Beta1 - 1).*cos(2.*phi) + (3 + Beta1).*...
                (B + A.*csc(theta).^2)));
            Term2 = ((C.^2).*((Beta2 - 1).^2).*I2.*(sin(2.*phi).^2))...
                ./(Beta2.*(B + A.*(csc(theta).^2)).*(C.*(Beta2 - 1).*cos(2.*phi) + (3 + Beta2).*...
                (B + A.*csc(theta).^2)));
            Term3 = ((C.^2).*((Beta3 - 1).^2).*I3.*(cos(2.*phi).^2))...
                ./(Beta3.*(B + A.*(csc(theta).^2)).*((3 + Beta3).*(B + A.*csc(theta).^2)...
                - C.*(Beta3 - 1).*sin(2.*phi)));
            Term4 = ((C.^2).*((Beta4 - 1).^2).*I4.*(cos(2.*phi).^2))...
                ./(Beta4.*(B + A.*(csc(theta).^2)).*((3 + Beta4).*(B + A.*csc(theta).^2)...
                + C.*(Beta4 - 1).*sin(2.*phi)));
            if n_detectors == 3
                Jval = time.*(Term1 + Term2 + Term3);
            elseif n_detectors == 4
                Jval = time.*(Term1 + Term2 + Term3 + Term4);
            end
            Jval = sum(Jval, 2);
        end
        
        function Jval = JphiphifNB(~, theta, phi, time, I, A, B, C, n_detectors)
            % time is time taken to observe photons
            % wavelength dependent), should be n_detector x wavelength
            % matrix
            % I is Intensity (can be wavelength dependent),
            % should be n_detector x wavelength matrix
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            
            if size(I, 2) == 1 % if intensity is a simple scalar
                if length(A) ~= 1 % if we have multiple wavelengths
                    I = I./length(A); % spread intensity across multiple wavelengths
                end
            end
            
            if n_detectors == 3
                I1 = I(1, :);
                I2 = I(2, :);
                I3 = I(3, :);
                I4 = 0;
            else
                I1 = I(1, :);
                I2 = I(2, :);
                I3 = I(3, :);
                I4 = I(4, :);
            end

            if n_detectors == 3
                Jval = time.*((1./(B + A.*(csc(theta).^2))).*(C.^2).*((sin(2.*phi).^2).*((I1./...
                    (B + C.*cos(2.*phi) + A.*(csc(theta).^2))) + (I3./...
                    (B - C.*cos(2.*phi) + A.*(csc(theta).^2)))) + (cos(2.*phi).^2).*...
                    (I2./(B + C.*sin(2.*phi) + A.*(csc(theta).^2)))));
            elseif n_detectors == 4
                Jval = time.*((1./(B + A.*(csc(theta).^2))).*(C.^2).*((sin(2.*phi).^2).*((I1./...
                    (B + C.*cos(2.*phi) + A.*(csc(theta).^2))) + (I3./...
                    (B - C.*cos(2.*phi) + A.*(csc(theta).^2)))) + (cos(2.*phi).^2).*...
                    (I2./(B + C.*sin(2.*phi) + A.*(csc(theta).^2))) + ...
                    (I4./(B - C.*sin(2.*phi) + A.*(csc(theta).^2)))));
            end
            Jval = sum(Jval, 2);
        end
       
        function NF = NormFactor(~, theta, a11, a13, a33, A, B, H)
            NF = (2.*(A + B).*(3 + H) - (A + 3.*A.*H + 4.*B.*(1 + H)).*(sin(theta).^2) + B.*(1 + 3.*H)...
                .*(sin(theta).^4)).*a11 + B.*(1 + 3.*H).*(cos(theta).^2).*(sin(theta).^2).*a13 ...
                - 0.25.*(2.*A + B - B.*cos(2.*theta)).*(-5 + H + (1 + 3.*H).*cos(2.*theta)).*a33;
        end
                
        function thetadepf = thetadepgetscatter(obj, a11, a13, a33, A, B, H)
            thetas = linspace(0, pi/2, 100000)';
            thetadepmean = mean(obj.NormFactor(thetas, a11, a13, a33, A, B, H), 2);
            m0 = trapz(sin(thetas), thetadepmean.*sin(thetas));
            m1 = trapz(sin(thetas), (thetadepmean.*sin(thetas)).*thetas)/m0;
            
            thetadep = thetadepmean./mean(obj.NormFactor(m1, a11, a13, a33, A, B, H), 2);   
            thetadepf = @(theta) interp1(thetas, thetadep, theta, 'v5cubic');
        end

        
        function NF = NormFactorF(~, theta, A, B)
            NF = 4.*(A + B.*(sin(theta).^2));
        end
        
%         function thetadepf = thetadepgetfluo(obj, A, B)
%             % 0.3183098861837907*pi is mean of sine distribution, see
%             % Mathematica notebook. However, shape of dependence on theta
%             % is related to a1, a13 and a3 in a complicated fashion
%             finalthets = zeros(1,10);
%             for i=1:10
%                 thetadepmean = @(thetval) mean(obj.NormFactorF(thetval, A, B), 2);
% 
%                 nsamples = 100000;
%                 angles = acos(rand(nsamples, 1)); % theta is from a sine distribution
% 
%                 fun = @(theta) abs(mean(mean(obj.NormFactorF(angles, A, B), 2))...
%                     - thetadepmean(theta));
%                 x0 = 1;
%                 options = optimoptions('fmincon','Display','off','TolFun',1e-9);
%                 finalthets(i) = fmincon(fun,x0,[],[],[],[],[],[],[],options);
%             end
%             finalthet = mean(finalthets);
%             x = linspace(0, pi/2, 10000)';
%             thetadep = mean(obj.NormFactorF(x, A, B), 2)...
%                 ./mean(obj.NormFactorF(finalthet, A, B), 2);
%             thetadepf = @(theta) interp1(x, thetadep, theta, 'v5cubic');
%         end
        
        function thetadepf = thetadepgetfluo(obj, A, B)
            thetas = linspace(0, pi/2, 100000)';
            thetadepmean = mean(obj.NormFactorF(thetas, A, B), 2);
            m0 = trapz(sin(thetas), thetadepmean.*sin(thetas));
            m1 = trapz(sin(thetas), (thetadepmean.*sin(thetas)).*thetas)/m0;
            
            thetadep = thetadepmean./mean(obj.NormFactorF(m1, A, B), 2);   
            thetadepf = @(theta) interp1(thetas, thetadep, theta, 'v5cubic');
        end
        
        function Qval = Qthetatheta1(~, theta, phi, Beta1, I1, A, B, C, H, a11, a13, a33)
            % time is time taken to observe photons
            % Beta is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent)
            % I is Intensity/(1-SBR^{-1}) (can be wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Qval = (C.^2).*(cos(2.*phi).^2).*((-2.*(16.*B.*((3 + H).^2) + A.*(147 + H.*(114 + 43.*H)) + ...
                4.*(1 + 3.*H).*(4.*B.*(3 + H) + A.*(11 + H)).*cos(2.*theta) + A.*((1 + 3.*H).^2).*cos(4.*theta))...
                .*sin(2.*theta).*(a11.^2) + 2.*A.*sin(2.*theta).*a33.*(-(1 + 3.*H).*(4.*(-5 + H).*cos(2.*theta) +...
                (1 + 3.*H).*(3 + cos(4.*theta))).*a13 + 2.*((-5 + H + (1 + 3.*H).*cos(2.*theta)).^2).*a33) + ...
                a11.*((1 + 3.*H).*(5.*A.*(1 + 3.*H).*sin(2.*theta) + 4.*(4.*B.*(3 + H) + A.*(11 + H)).*sin(4.*theta) +...
                A.*(1 + 3.*H).*sin(6.*theta)).*a13 + 16.*(2.*(3 + H).*(2.*A.*(1 + H) + B.*(3 + H)).*sin(2.*theta) ...
                - (1 + 3.*H).*(4.*A + B.*(3 + H)).*sin(4.*theta)).*a33)).^2).*(((1./Beta1) - 1).^2).*I1;
        end
        
        function Qval = Qthetatheta1NB(~, theta, phi, I1, A, B, C, H, a11, a13, a33)
            % time is time taken to observe photons
            % I is Intensity (can be wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Qval = (C.^2).*(cos(2.*phi).^2).*((-2.*(16.*B.*((3 + H).^2) + A.*(147 + H.*(114 + 43.*H)) + ...
                4.*(1 + 3.*H).*(4.*B.*(3 + H) + A.*(11 + H)).*cos(2.*theta) + A.*((1 + 3.*H).^2).*cos(4.*theta))...
                .*sin(2.*theta).*(a11.^2) + 2.*A.*sin(2.*theta).*a33.*(-(1 + 3.*H).*(4.*(-5 + H).*cos(2.*theta) +...
                (1 + 3.*H).*(3 + cos(4.*theta))).*a13 + 2.*((-5 + H + (1 + 3.*H).*cos(2.*theta)).^2).*a33) + ...
                a11.*((1 + 3.*H).*(5.*A.*(1 + 3.*H).*sin(2.*theta) + 4.*(4.*B.*(3 + H) + A.*(11 + H)).*sin(4.*theta) +...
                A.*(1 + 3.*H).*sin(6.*theta)).*a13 + 16.*(2.*(3 + H).*(2.*A.*(1 + H) + B.*(3 + H)).*sin(2.*theta) ...
                - (1 + 3.*H).*(4.*A + B.*(3 + H)).*sin(4.*theta)).*a33)).^2).*I1;
        end
        
        function Mval = Mthetatheta1(~, theta, phi, Beta1, A, B, C, H, NF, a11, a13, a33)
            % time is time taken to observe photons
            % Beta is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Mval = 4096.*(NF.^4).*((1./(16.*NF)).*(4.*(2.*(A + B).*(3 + H) - (A + 3.*A.*H + 4.*B.*(1 + H) + ...
                4.*(1 + H).*C.*cos(2.*phi)).*(sin(theta).^2) + (1 + 3.*H).*(B + C.*cos(2.*phi)).*(sin(theta).^4)).*a11...
                + (1 + 3.*H).*(B + C.*cos(2.*phi)).*(sin(2.*theta).^2).*a13 + 4.*(2 - 2.*H + (1 + 3.*H).*(sin(theta).^2))...
                .*(A + (B + C.*cos(2.*phi)).*(sin(theta).^2)).*a33).*(1 - (1./Beta1)) + (1./Beta1));
        end
        
        function Mval = Mthetatheta1NB(~, theta, phi, A, B, C, H, NF, a11, a13, a33)
            % time is time taken to observe photons
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Mval = 4096.*(NF.^4).*((1./(16.*NF)).*(4.*(2.*(A + B).*(3 + H) - (A + 3.*A.*H + 4.*B.*(1 + H) + ...
                4.*(1 + H).*C.*cos(2.*phi)).*(sin(theta).^2) + (1 + 3.*H).*(B + C.*cos(2.*phi)).*(sin(theta).^4)).*a11...
                + (1 + 3.*H).*(B + C.*cos(2.*phi)).*(sin(2.*theta).^2).*a13 + 4.*(2 - 2.*H + (1 + 3.*H).*(sin(theta).^2))...
                .*(A + (B + C.*cos(2.*phi)).*(sin(theta).^2)).*a33));
        end

        function Qval = Qthetatheta2(~, theta, phi, Beta2, I2, A, B, C, H, a11, a13, a33)
            % Beta_n is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent)
            % I_n is Intensity/(1-SBR^{-1}) (can be wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms

            Qval = (C.^2).*(sin(2.*theta).^2).*(sin(2.*phi).^2).*((-(16.*B.*((3 + H).^2) + A.*(147 + H.*(114 + 43.*H))...
                + 4.*(1 + 3.*H).*(4.*B.*(3 + H) + A.*(11 + H)).*cos(2.*theta) + A.*((1 + 3.*H).^2).*cos(4.*theta))...
                .*(a11.^2) + A.*a33.*(-(1 + 3.*H).*(4.*(-5 + H).*cos(2.*theta) + (1 + 3.*H).*(3 + cos(4.*theta))).*a13...
                + 2.*((-5 + H + (1 + 3.*H).*cos(2.*theta)).^2).*a33) + a11.*((1 + 3.*H).*(4.*(4.*B.*(3 + H) + A.*(11 + H))...
                .*cos(2.*theta) + A.*(1 + 3.*H).*(3 + cos(4.*theta))).*a13 + 16.*((3 + H).*(2.*A.*(1 + H) + B.*(3 + H))...
                - (1 + 3.*H).*(4.*A + B.*(3 + H)).*cos(2.*theta)).*a33)).^2).*((((1./Beta2) - 1)).^2).*I2;
        end
        
        function Qval = Qthetatheta2NB(~, theta, phi, I2, A, B, C, H, a11, a13, a33)
            % I_n is Intensity/(1-SBR^{-1}) (can be wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms

            Qval = (C.^2).*(sin(2.*theta).^2).*(sin(2.*phi).^2).*((-(16.*B.*((3 + H).^2) + A.*(147 + H.*(114 + 43.*H))...
                + 4.*(1 + 3.*H).*(4.*B.*(3 + H) + A.*(11 + H)).*cos(2.*theta) + A.*((1 + 3.*H).^2).*cos(4.*theta))...
                .*(a11.^2) + A.*a33.*(-(1 + 3.*H).*(4.*(-5 + H).*cos(2.*theta) + (1 + 3.*H).*(3 + cos(4.*theta))).*a13...
                + 2.*((-5 + H + (1 + 3.*H).*cos(2.*theta)).^2).*a33) + a11.*((1 + 3.*H).*(4.*(4.*B.*(3 + H) + A.*(11 + H))...
                .*cos(2.*theta) + A.*(1 + 3.*H).*(3 + cos(4.*theta))).*a13 + 16.*((3 + H).*(2.*A.*(1 + H) + B.*(3 + H))...
                - (1 + 3.*H).*(4.*A + B.*(3 + H)).*cos(2.*theta)).*a33)).^2).*I2;
        end
        
        function Mval = Mthetatheta2(~, theta, phi, Beta2, A, B, C, H, NF, a11, a13, a33)
            % time is time taken to observe photons
            % Beta is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent), should be n_detector x wavelength
            % matrix
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Mval = 1024.*(NF.^4).*((1./(16.*NF)).*(4.*(2.*(A + B).*(3 + H) + (1 + 3.*H).*(sin(theta).^4).*...
                (B + C.*sin(2.*phi)) - (sin(theta).^2).*(A + 3.*A.*H + 4.*B.*(1 + H) + 4.*(1 + H).*C.*sin(2.*phi))).*a11...
                + (1 + 3.*H).*(sin(2.*theta).^2).*(B + C.*sin(2.*phi)).*a13 + 4.*(2 - 2.*H + (1 + 3.*H).*(sin(theta).^2))...
                .*(A + (sin(theta).^2).*(B + C.*sin(2.*phi))).*a33).*(1 - (1./Beta2)) + (1./Beta2));
        end
        
        function Mval = Mthetatheta2NB(~, theta, phi, A, B, C, H, NF, a11, a13, a33)
            % time is time taken to observe photons
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Mval = 1024.*(NF.^4).*((1./(16.*NF)).*(4.*(2.*(A + B).*(3 + H) + (1 + 3.*H).*(sin(theta).^4).*...
                (B + C.*sin(2.*phi)) - (sin(theta).^2).*(A + 3.*A.*H + 4.*B.*(1 + H) + 4.*(1 + H).*C.*sin(2.*phi))).*a11...
                + (1 + 3.*H).*(sin(2.*theta).^2).*(B + C.*sin(2.*phi)).*a13 + 4.*(2 - 2.*H + (1 + 3.*H).*(sin(theta).^2))...
                .*(A + (sin(theta).^2).*(B + C.*sin(2.*phi))).*a33));
        end

        function Qval = Qthetatheta3(~, theta, phi, Beta3, I3, A, B, C, H, a11, a13, a33)
            % Beta_n is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent)
            % I_n is Intensity/(1-SBR^{-1}) (can be wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Qval = (C.^2).*(cos(2.*phi).^2).*((2.*(16.*B.*((3 + H).^2) + A.*(147 + H.*(114 + 43.*H)) + ...
                4.*(1 + 3.*H).*(4.*B.*(3 + H) + A.*(11 + H)).*cos(2.*theta) + A.*((1 + 3.*H).^2).*cos(4.*theta))...
                .*sin(2.*theta).*(a11.^2) + 2.*A.*sin(2.*theta).*a33.*((1 + 3.*H).*(4.*(-5 + H).*cos(2.*theta) +...
                (1 + 3.*H).*(3 + cos(4.*theta))).*a13 - 2.*((-5 + H + (1 + 3.*H).*cos(2.*theta)).^2).*a33) + ...
                a11.*(-2.*(1 + 3.*H).*(4.*(4.*B.*(3 + H) + A.*(11 + H)).*cos(2.*theta) + A.*(1 + 3.*H)...
                .*(3 + cos(4.*theta))).*sin(2.*theta).*a13 + 16.*(-2.*(3 + H).*(2.*A.*(1 + H) + B.*(3 + H))...
                .*sin(2.*theta) + (1 + 3.*H).*(4.*A + B.*(3 + H)).*sin(4.*theta)).*a33)).^2).*(((1./Beta3) - 1).^2).*I3;
        end
        
        function Qval = Qthetatheta3NB(~, theta, phi, I3, A, B, C, H, a11, a13, a33)
            % I_n is Intensity/(1-SBR^{-1}) (can be wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Qval = (C.^2).*(cos(2.*phi).^2).*((2.*(16.*B.*((3 + H).^2) + A.*(147 + H.*(114 + 43.*H)) + ...
                4.*(1 + 3.*H).*(4.*B.*(3 + H) + A.*(11 + H)).*cos(2.*theta) + A.*((1 + 3.*H).^2).*cos(4.*theta))...
                .*sin(2.*theta).*(a11.^2) + 2.*A.*sin(2.*theta).*a33.*((1 + 3.*H).*(4.*(-5 + H).*cos(2.*theta) +...
                (1 + 3.*H).*(3 + cos(4.*theta))).*a13 - 2.*((-5 + H + (1 + 3.*H).*cos(2.*theta)).^2).*a33) + ...
                a11.*(-2.*(1 + 3.*H).*(4.*(4.*B.*(3 + H) + A.*(11 + H)).*cos(2.*theta) + A.*(1 + 3.*H)...
                .*(3 + cos(4.*theta))).*sin(2.*theta).*a13 + 16.*(-2.*(3 + H).*(2.*A.*(1 + H) + B.*(3 + H))...
                .*sin(2.*theta) + (1 + 3.*H).*(4.*A + B.*(3 + H)).*sin(4.*theta)).*a33)).^2).*I3;
        end

        function Mval = Mthetatheta3(~, theta, phi, Beta3, A, B, C, H, NF, a11, a13, a33)
            % time is time taken to observe photons
            % Beta is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Mval = 4096.*(NF.^4).*((1./(16.*NF)).*(4.*(2.*(A + B).*(3 + H) - (A + 3.*A.*H + 4.*B.*(1 + H) - ...
                4.*(1 + H).*C.*cos(2.*phi)).*(sin(theta).^2) + (1 + 3.*H).*(B - C.*cos(2.*phi)).*(sin(theta).^4)).*a11...
                + (1 + 3.*H).*(B - C.*cos(2.*phi)).*(sin(2.*theta).^2).*a13 + 4.*(2 - 2.*H + (1 + 3.*H).*(sin(theta).^2))...
                .*(A + (B - C.*cos(2.*phi)).*(sin(theta).^2)).*a33).*(1 - (1./Beta3)) + (1./Beta3));
        end

        function Mval = Mthetatheta3NB(~, theta, phi, A, B, C, H, NF, a11, a13, a33)
            % time is time taken to observe photons
            % Beta is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Mval = 4096.*(NF.^4).*((1./(16.*NF)).*(4.*(2.*(A + B).*(3 + H) - (A + 3.*A.*H + 4.*B.*(1 + H) - ...
                4.*(1 + H).*C.*cos(2.*phi)).*(sin(theta).^2) + (1 + 3.*H).*(B - C.*cos(2.*phi)).*(sin(theta).^4)).*a11...
                + (1 + 3.*H).*(B - C.*cos(2.*phi)).*(sin(2.*theta).^2).*a13 + 4.*(2 - 2.*H + (1 + 3.*H).*(sin(theta).^2))...
                .*(A + (B - C.*cos(2.*phi)).*(sin(theta).^2)).*a33));
        end
        
        function Qval = Qthetatheta4(~, theta, phi, Beta4, I4, A, B, C, H, a11, a13, a33)
            % Beta_n is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent)
            % I_n is Intensity/(1-SBR^{-1}) (can be wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Qval = (C.^2).*(sin(2.*theta).^2).*(sin(2.*phi).^2).*(((16.*B.*((3 + H).^2) + A.*(147 + H.*(114 + 43.*H))...
                + (1 + 3.*H).*(4.*(4.*B.*(3 + H) + A.*(11 + H)).*cos(2.*theta) + A.*(1 + 3.*H).*cos(4.*theta)))...
                .*(a11.^2) + A.*a33.*((1 + 3.*H).*(4.*(-5 + H).*cos(2.*theta) + (1 + 3.*H).*(3 + cos(4.*theta))).*a13...
                - 2.*((-5 + H + (1 + 3.*H).*cos(2.*theta)).^2).*a33) + a11.*(-(1 + 3.*H).*(4.*(4.*B.*(3 + H) + A.*(11 + H))...
                .*cos(2.*theta) + A.*(1 + 3.*H).*(3 + cos(4.*theta))).*a13 + 16.*(-(3 + H).*(2.*A.*(1 + H) + B.*(3 + H))...
                + (1 + 3.*H).*(4.*A + B.*(3 + H)).*cos(2.*theta)).*a33)).^2).*((((1./Beta4) - 1)).^2).*I4;
        end
        
        function Qval = Qthetatheta4NB(~, theta, phi, I4, A, B, C, H, a11, a13, a33)
            % I_n is Intensity/(1-SBR^{-1}) (can be wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Qval = (C.^2).*(sin(2.*theta).^2).*(sin(2.*phi).^2).*(((16.*B.*((3 + H).^2) + A.*(147 + H.*(114 + 43.*H))...
                + (1 + 3.*H).*(4.*(4.*B.*(3 + H) + A.*(11 + H)).*cos(2.*theta) + A.*(1 + 3.*H).*cos(4.*theta)))...
                .*(a11.^2) + A.*a33.*((1 + 3.*H).*(4.*(-5 + H).*cos(2.*theta) + (1 + 3.*H).*(3 + cos(4.*theta))).*a13...
                - 2.*((-5 + H + (1 + 3.*H).*cos(2.*theta)).^2).*a33) + a11.*(-(1 + 3.*H).*(4.*(4.*B.*(3 + H) + A.*(11 + H))...
                .*cos(2.*theta) + A.*(1 + 3.*H).*(3 + cos(4.*theta))).*a13 + 16.*(-(3 + H).*(2.*A.*(1 + H) + B.*(3 + H))...
                + (1 + 3.*H).*(4.*A + B.*(3 + H)).*cos(2.*theta)).*a33)).^2).*I4;
        end

        function Mval = Mthetatheta4(~, theta, phi, Beta4, A, B, C, H, NF, a11, a13, a33)
            % time is time taken to observe photons
            % Beta is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent)
            % I is Intensity/(1-SBR^{-1}) (can be wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Mval = 1024.*(NF.^4).*((1./(16.*NF)).*(4.*(2.*(A + B).*(3 + H) + (1 + 3.*H).*(sin(theta).^4).*...
                (B - C.*sin(2.*phi)) - (sin(theta).^2).*(A + 3.*A.*H + 4.*B.*(1 + H) - 4.*(1 + H).*C.*sin(2.*phi))).*a11...
                + (1 + 3.*H).*(sin(2.*theta).^2).*(B - C.*sin(2.*phi)).*a13 + 4.*(2 - 2.*H + (1 + 3.*H).*(sin(theta).^2))...
                .*(A + (sin(theta).^2).*(B - C.*sin(2.*phi))).*a33).*(1 - (1./Beta4)) + (1./Beta4));
        end
        
        function Mval = Mthetatheta4NB(~, theta, phi, A, B, C, H, NF, a11, a13, a33)
            % I is Intensity (can be wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Mval = 1024.*(NF.^4).*((1./(16.*NF)).*(4.*(2.*(A + B).*(3 + H) + (1 + 3.*H).*(sin(theta).^4).*...
                (B - C.*sin(2.*phi)) - (sin(theta).^2).*(A + 3.*A.*H + 4.*B.*(1 + H) - 4.*(1 + H).*C.*sin(2.*phi))).*a11...
                + (1 + 3.*H).*(sin(2.*theta).^2).*(B - C.*sin(2.*phi)).*a13 + 4.*(2 - 2.*H + (1 + 3.*H).*(sin(theta).^2))...
                .*(A + (sin(theta).^2).*(B - C.*sin(2.*phi))).*a33));
        end

        function Qval = Qthetaphi1(~, theta, phi, Beta1, I1, A, B, C, H, a11, a13, a33)
            % Beta_n is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent)
            % I_n is Intensity/(1-SBR^{-1}) (can be wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Qval = (C.^2).*cos(2.*phi).*(sin(theta).^2).*sin(2.*phi).*((7 + 5.*H + (1 + 3.*H).*cos(2.*theta)).*a11...
                - 2.*(1 + 3.*H).*(cos(theta).^2).*a13 + (-5 + H + (1 + 3.*H).*cos(2.*theta)).*a33).*...
                (-2.*(16.*B.*((3 + H).^2) + A.*(147 + H.*(114 + 43.*H)) + 4.*(1 + 3.*H).*(4.*B.*(3 + H) + ...
                A.*(11 + H)).*cos(2.*theta) + A.*((1 + 3.*H).^2).*cos(4.*theta)).*sin(2.*theta).*(a11.^2)...
                + 2.*A.*sin(2.*theta).*a33.*(-(1 + 3.*H).*(4.*(-5 + H).*cos(2.*theta) + (1 + 3.*H)...
                .*(3 + cos(4.*theta))).*a13 + 2.*((-5 + H + (1 + 3.*H).*cos(2.*theta)).^2).*a33) + ...
                a11.*((1 + 3.*H).*(5.*A.*(1 + 3.*H).*sin(2.*theta) + 4.*(4.*B.*(3 + H) + A.*(11 + H)).*...
                sin(4.*theta) + A.*(1 + 3.*H).*sin(6.*theta)).*a13 + ...
                16.*(2.*(3 + H).*(2.*A.*(1 + H) + B.*(3 + H)).*sin(2.*theta) - (1 + 3.*H).*(4.*A + B.*(3 + H))...
                .*sin(4.*theta)).*a33)).*(((1./Beta1) - 1).^2).*I1;
        end
        
        function Qval = Qthetaphi1NB(~, theta, phi, I1, A, B, C, H, a11, a13, a33)
            % I_n is Intensity (can be wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Qval = (C.^2).*cos(2.*phi).*(sin(theta).^2).*sin(2.*phi).*((7 + 5.*H + (1 + 3.*H).*cos(2.*theta)).*a11...
                - 2.*(1 + 3.*H).*(cos(theta).^2).*a13 + (-5 + H + (1 + 3.*H).*cos(2.*theta)).*a33).*...
                (-2.*(16.*B.*((3 + H).^2) + A.*(147 + H.*(114 + 43.*H)) + 4.*(1 + 3.*H).*(4.*B.*(3 + H) + ...
                A.*(11 + H)).*cos(2.*theta) + A.*((1 + 3.*H).^2).*cos(4.*theta)).*sin(2.*theta).*(a11.^2)...
                + 2.*A.*sin(2.*theta).*a33.*(-(1 + 3.*H).*(4.*(-5 + H).*cos(2.*theta) + (1 + 3.*H)...
                .*(3 + cos(4.*theta))).*a13 + 2.*((-5 + H + (1 + 3.*H).*cos(2.*theta)).^2).*a33) + ...
                a11.*((1 + 3.*H).*(5.*A.*(1 + 3.*H).*sin(2.*theta) + 4.*(4.*B.*(3 + H) + A.*(11 + H)).*...
                sin(4.*theta) + A.*(1 + 3.*H).*sin(6.*theta)).*a13 + ...
                16.*(2.*(3 + H).*(2.*A.*(1 + H) + B.*(3 + H)).*sin(2.*theta) - (1 + 3.*H).*(4.*A + B.*(3 + H))...
                .*sin(4.*theta)).*a33)).*I1;
        end

        function Mval = Mthetaphi1(~, theta, phi, Beta1, A, B, C, H, NF, a11, a13, a33)
            % time is time taken to observe photons
            % Beta_n is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent)
            % I_n is Intensity/(1-SBR^{-1}) (can be wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Mval = 16.*(NF.^2).*((4.*(2.*(A + B).*(3 + H) - (A + 3.*A.*H + 4.*B.*(1 + H) + ...
                4.*(1 + H).*C.*cos(2.*phi)).*(sin(theta).^2) + (1 + 3.*H).*(B + C.*cos(2.*phi)).*(sin(theta).^4)).*a11...
                + (1 + 3.*H).*(B + C.*cos(2.*phi)).*(sin(2.*theta).^2).*a13 + ...
                4.*(2 - 2.*H + (1 + 3.*H).*(sin(theta).^2)).*(A + (B + C.*cos(2.*phi)).*(sin(theta).^2)).*a33)...
                .*(1 - (1./Beta1)) + ((16.*NF)./Beta1));
        end
        
        function Mval = Mthetaphi1NB(~, theta, phi, A, B, C, H, NF, a11, a13, a33)
            % I_n is Intensity (can be wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Mval = 16.*(NF.^2).*((4.*(2.*(A + B).*(3 + H) - (A + 3.*A.*H + 4.*B.*(1 + H) + ...
                4.*(1 + H).*C.*cos(2.*phi)).*(sin(theta).^2) + (1 + 3.*H).*(B + C.*cos(2.*phi)).*(sin(theta).^4)).*a11...
                + (1 + 3.*H).*(B + C.*cos(2.*phi)).*(sin(2.*theta).^2).*a13 + ...
                4.*(2 - 2.*H + (1 + 3.*H).*(sin(theta).^2)).*(A + (B + C.*cos(2.*phi)).*(sin(theta).^2)).*a33));
        end
        
        function Qval = Qthetaphi2(~, theta, phi, Beta2, I2, A, B, C, H, a11, a13, a33)
            % Beta_n is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent)
            % I_n is Intensity/(1-SBR^{-1}) (can be wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Qval = (C.^2).*cos(2.*phi).*(sin(theta).^2).*sin(2.*theta).*sin(2.*phi).*((7 + 5.*H + (1 + 3.*H)...
                .*cos(2.*theta)).*a11 - 2.*(1 + 3.*H).*(cos(theta).^2).*a13 + (-5 + H + (1 + 3.*H).*cos(2.*theta))...
                .*a33).*(-(16.*B.*((3 + H).^2) + A.*(147 + H.*(114 + 43.*H)) + 4.*(1 + 3.*H).*(4.*B.*(3 + H)...
                + A.*(11 + H)).*cos(2.*theta) + A.*((1 + 3.*H).^2).*cos(4.*theta)).*(a11.^2) + A.*a33.*...
                (-(1 + 3.*H).*(4.*(-5 + H).*cos(2.*theta) + (1 + 3.*H).*(3 + cos(4.*theta))).*a13 + ...
                2.*((-5 + H + (1 + 3.*H).*cos(2.*theta)).^2).*a33) + ...
                a11.*((1 + 3.*H).*(4.*(4.*B.*(3 + H) + A.*(11 + H)).*cos(2.*theta) + A.*(1 + 3.*H).*...
                (3 + cos(4.*theta))).*a13 + 16.*((3 + H).*(2.*A.*(1 + H) + B.*(3 + H)) ...
                - (1 + 3.*H).*(4.*A + B.*(3 + H)).*cos(2.*theta)).*a33)).*(((1./Beta2) - 1).^2).*I2;
        end
        
        function Qval = Qthetaphi2NB(~, theta, phi, I2, A, B, C, H, a11, a13, a33)
            % I_n is Intensity (can be wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Qval = (C.^2).*cos(2.*phi).*(sin(theta).^2).*sin(2.*theta).*sin(2.*phi).*((7 + 5.*H + (1 + 3.*H)...
                .*cos(2.*theta)).*a11 - 2.*(1 + 3.*H).*(cos(theta).^2).*a13 + (-5 + H + (1 + 3.*H).*cos(2.*theta))...
                .*a33).*(-(16.*B.*((3 + H).^2) + A.*(147 + H.*(114 + 43.*H)) + 4.*(1 + 3.*H).*(4.*B.*(3 + H)...
                + A.*(11 + H)).*cos(2.*theta) + A.*((1 + 3.*H).^2).*cos(4.*theta)).*(a11.^2) + A.*a33.*...
                (-(1 + 3.*H).*(4.*(-5 + H).*cos(2.*theta) + (1 + 3.*H).*(3 + cos(4.*theta))).*a13 + ...
                2.*((-5 + H + (1 + 3.*H).*cos(2.*theta)).^2).*a33) + ...
                a11.*((1 + 3.*H).*(4.*(4.*B.*(3 + H) + A.*(11 + H)).*cos(2.*theta) + A.*(1 + 3.*H).*...
                (3 + cos(4.*theta))).*a13 + 16.*((3 + H).*(2.*A.*(1 + H) + B.*(3 + H)) ...
                - (1 + 3.*H).*(4.*A + B.*(3 + H)).*cos(2.*theta)).*a33)).*I2;
        end

        
        function Mval = Mthetaphi2(~, theta, phi, Beta2, A, B, C, H, NF, a11, a13, a33)
            % Beta_n is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent)
            % I_n is Intensity/(1-SBR^{-1}) (can be wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Mval = 8.*(NF.^2).*((4.*(2.*(A + B).*(3 + H) + (1 + 3.*H).*(sin(theta).^4).*(B + C.*sin(2.*phi))...
                - (sin(theta).^2).*(A + 3.*A.*H + 4.*B.*(1 + H) + 4.*(1 + H).*C.*sin(2.*phi))).*a11 ...
                + (1 + 3.*H).*(sin(2.*theta).^2).*(B + C.*sin(2.*phi)).*a13 + 4.*(2 - 2.*H + (1 + 3.*H)...
                .*(sin(theta).^2)).*(A + (sin(theta).^2).*(B + C.*sin(2.*phi))).*a33).*(1 - (1./Beta2)) + ...
                ((16.*NF)./Beta2));
        end
        
        function Mval = Mthetaphi2NB(~, theta, phi, A, B, C, H, NF, a11, a13, a33)
            % I_n is Intensity (can be wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Mval = 8.*(NF.^2).*((4.*(2.*(A + B).*(3 + H) + (1 + 3.*H).*(sin(theta).^4).*(B + C.*sin(2.*phi))...
                - (sin(theta).^2).*(A + 3.*A.*H + 4.*B.*(1 + H) + 4.*(1 + H).*C.*sin(2.*phi))).*a11 ...
                + (1 + 3.*H).*(sin(2.*theta).^2).*(B + C.*sin(2.*phi)).*a13 + 4.*(2 - 2.*H + (1 + 3.*H)...
                .*(sin(theta).^2)).*(A + (sin(theta).^2).*(B + C.*sin(2.*phi))).*a33));
        end
        
        function Qval = Qthetaphi3(~, theta, phi, Beta3, I3, A, B, C, H, a11, a13, a33)
            % Beta_n is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent)
            % I_n is Intensity/(1-SBR^{-1}) (can be wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Qval = (C.^2).*cos(2.*phi).*(sin(theta).^2).*sin(2.*phi).*((7 + 5.*H + (1 + 3.*H).*cos(2.*theta)).*a11...
                - 2.*(1 + 3.*H).*(cos(theta).^2).*a13 + (-5 + H + (1 + 3.*H).*cos(2.*theta)).*a33).*...
                (2.*(16.*B.*((3 + H).^2) + A.*(147 + H.*(114 + 43.*H)) + 4.*(1 + 3.*H).*(4.*B.*(3 + H) + ...
                A.*(11 + H)).*cos(2.*theta) + A.*((1 + 3.*H).^2).*cos(4.*theta)).*sin(2.*theta).*(a11.^2)...
                + 2.*A.*sin(2.*theta).*a33.*((1 + 3.*H).*(4.*(-5 + H).*cos(2.*theta) + (1 + 3.*H)...
                .*(3 + cos(4.*theta))).*a13 - 2.*((-5 + H + (1 + 3.*H).*cos(2.*theta)).^2).*a33) + ...
                a11.*(-2.*(1 + 3.*H).*(4.*(4.*B.*(3 + H) + A.*(11 + H)).*cos(2.*theta) + A.*(1 + 3.*H)...
                .*(3 + cos(4.*theta))).*sin(2.*theta).*a13 + 16.*(-2.*(3 + H).*(2.*A.*(1 + H)...
                + B.*(3 + H)).*sin(2.*theta) + (1 + 3.*H).*(4.*A + B.*(3 + H)).*sin(4.*theta)).*a33))...
                .*(((1./Beta3) - 1).^2).*I3;
        end
        
        function Qval = Qthetaphi3NB(~, theta, phi, I3, A, B, C, H, a11, a13, a33)
            % I_n is Intensity (can be wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Qval = (C.^2).*cos(2.*phi).*(sin(theta).^2).*sin(2.*phi).*((7 + 5.*H + (1 + 3.*H).*cos(2.*theta)).*a11...
                - 2.*(1 + 3.*H).*(cos(theta).^2).*a13 + (-5 + H + (1 + 3.*H).*cos(2.*theta)).*a33).*...
                (2.*(16.*B.*((3 + H).^2) + A.*(147 + H.*(114 + 43.*H)) + 4.*(1 + 3.*H).*(4.*B.*(3 + H) + ...
                A.*(11 + H)).*cos(2.*theta) + A.*((1 + 3.*H).^2).*cos(4.*theta)).*sin(2.*theta).*(a11.^2)...
                + 2.*A.*sin(2.*theta).*a33.*((1 + 3.*H).*(4.*(-5 + H).*cos(2.*theta) + (1 + 3.*H)...
                .*(3 + cos(4.*theta))).*a13 - 2.*((-5 + H + (1 + 3.*H).*cos(2.*theta)).^2).*a33) + ...
                a11.*(-2.*(1 + 3.*H).*(4.*(4.*B.*(3 + H) + A.*(11 + H)).*cos(2.*theta) + A.*(1 + 3.*H)...
                .*(3 + cos(4.*theta))).*sin(2.*theta).*a13 + 16.*(-2.*(3 + H).*(2.*A.*(1 + H)...
                + B.*(3 + H)).*sin(2.*theta) + (1 + 3.*H).*(4.*A + B.*(3 + H)).*sin(4.*theta)).*a33))...
                .*I3;
        end

        function Mval = Mthetaphi3(~, theta, phi, Beta3, A, B, C, H, NF, a11, a13, a33)
            % Beta_n is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent)
            % I_n is Intensity/(1-SBR^{-1}) (can be wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Mval = 16.*(NF.^2).*((4.*(2.*(A + B).*(3 + H) - (A + 3.*A.*H + 4.*B.*(1 + H) - ...
                4.*(1 + H).*C.*cos(2.*phi)).*(sin(theta).^2) + (1 + 3.*H).*(B - C.*cos(2.*phi)).*(sin(theta).^4)).*a11...
                + (1 + 3.*H).*(B - C.*cos(2.*phi)).*(sin(2.*theta).^2).*a13 + ...
                4.*(2 - 2.*H + (1 + 3.*H).*(sin(theta).^2)).*(A + (B - C.*cos(2.*phi)).*(sin(theta).^2)).*a33)...
                .*(1 - (1./Beta3)) + ((16.*NF)./Beta3));
        end
        
        function Mval = Mthetaphi3NB(~, theta, phi, A, B, C, H, NF, a11, a13, a33)
            % I_n is Intensity (can be wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Mval = 16.*(NF.^2).*((4.*(2.*(A + B).*(3 + H) - (A + 3.*A.*H + 4.*B.*(1 + H) - ...
                4.*(1 + H).*C.*cos(2.*phi)).*(sin(theta).^2) + (1 + 3.*H).*(B - C.*cos(2.*phi)).*(sin(theta).^4)).*a11...
                + (1 + 3.*H).*(B - C.*cos(2.*phi)).*(sin(2.*theta).^2).*a13 + ...
                4.*(2 - 2.*H + (1 + 3.*H).*(sin(theta).^2)).*(A + (B - C.*cos(2.*phi)).*(sin(theta).^2)).*a33));
        end
        
        function Qval = Qthetaphi4(~, theta, phi, Beta4, I4, A, B, C, H, a11, a13, a33)
            % Beta_n is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent)
            % I_n is Intensity/(1-SBR^{-1}) (can be wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Qval = (C.^2).*cos(2.*phi).*(sin(theta).^2).*sin(2.*theta).*sin(2.*phi).*((7 + 5.*H + (1 + 3.*H)...
                .*cos(2.*theta)).*a11 - 2.*(1 + 3.*H).*(cos(theta).^2).*a13 + (-5 + H + (1 + 3.*H).*cos(2.*theta))...
                .*a33).*((16.*B.*((3 + H).^2) + A.*(147 + H.*(114 + 43.*H)) + (1 + 3.*H).*(4.*(4.*B.*(3 + H)...
                + A.*(11 + H)).*cos(2.*theta) + A.*(1 + 3.*H).*cos(4.*theta))).*(a11.^2) + A.*a33.*...
                ((1 + 3.*H).*(4.*(-5 + H).*cos(2.*theta) + (1 + 3.*H).*(3 + cos(4.*theta))).*a13 - ...
                2.*((-5 + H + (1 + 3.*H).*cos(2.*theta)).^2).*a33) + ...
                a11.*(-(1 + 3.*H).*(4.*(4.*B.*(3 + H) + A.*(11 + H)).*cos(2.*theta) + A.*(1 + 3.*H).*...
                (3 + cos(4.*theta))).*a13 + 16.*(-(3 + H).*(2.*A.*(1 + H) + B.*(3 + H)) ...
                + (1 + 3.*H).*(4.*A + B.*(3 + H)).*cos(2.*theta)).*a33)).*(((1./Beta4) - 1).^2).*I4;

        end
        
        function Qval = Qthetaphi4NB(~, theta, phi, I4, A, B, C, H, a11, a13, a33)
            % I_n is Intensity (can be wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Qval = (C.^2).*cos(2.*phi).*(sin(theta).^2).*sin(2.*theta).*sin(2.*phi).*((7 + 5.*H + (1 + 3.*H)...
                .*cos(2.*theta)).*a11 - 2.*(1 + 3.*H).*(cos(theta).^2).*a13 + (-5 + H + (1 + 3.*H).*cos(2.*theta))...
                .*a33).*((16.*B.*((3 + H).^2) + A.*(147 + H.*(114 + 43.*H)) + (1 + 3.*H).*(4.*(4.*B.*(3 + H)...
                + A.*(11 + H)).*cos(2.*theta) + A.*(1 + 3.*H).*cos(4.*theta))).*(a11.^2) + A.*a33.*...
                ((1 + 3.*H).*(4.*(-5 + H).*cos(2.*theta) + (1 + 3.*H).*(3 + cos(4.*theta))).*a13 - ...
                2.*((-5 + H + (1 + 3.*H).*cos(2.*theta)).^2).*a33) + ...
                a11.*(-(1 + 3.*H).*(4.*(4.*B.*(3 + H) + A.*(11 + H)).*cos(2.*theta) + A.*(1 + 3.*H).*...
                (3 + cos(4.*theta))).*a13 + 16.*(-(3 + H).*(2.*A.*(1 + H) + B.*(3 + H)) ...
                + (1 + 3.*H).*(4.*A + B.*(3 + H)).*cos(2.*theta)).*a33)).*I4;

        end
        
        function Mval = Mthetaphi4(~, theta, phi, Beta4, A, B, C, H, NF, a11, a13, a33)
            % time is time taken to observe photons
            % Beta_n is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent)
            % I_n is Intensity/(1-SBR^{-1}) (can be wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Mval = 8.*(NF.^2).*((4.*(2.*(A + B).*(3 + H) + (1 + 3.*H).*(sin(theta).^4).*(B - C.*sin(2.*phi))...
                - (sin(theta).^2).*(A + 3.*A.*H + 4.*B.*(1 + H) - 4.*(1 + H).*C.*sin(2.*phi))).*a11 ...
                + (1 + 3.*H).*(sin(2.*theta).^2).*(B - C.*sin(2.*phi)).*a13 + 4.*(2 - 2.*H + (1 + 3.*H)...
                .*(sin(theta).^2)).*(A + (sin(theta).^2).*(B - C.*sin(2.*phi))).*a33).*(1 - (1./Beta4)) + ...
                ((16.*NF)./Beta4));
        end
        
        function Mval = Mthetaphi4NB(~, theta, phi, A, B, C, H, NF, a11, a13, a33)
            % time is time taken to observe photons
            % Beta_n is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent)
            % I_n is Intensity/(1-SBR^{-1}) (can be wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Mval = 8.*(NF.^2).*((4.*(2.*(A + B).*(3 + H) + (1 + 3.*H).*(sin(theta).^4).*(B - C.*sin(2.*phi))...
                - (sin(theta).^2).*(A + 3.*A.*H + 4.*B.*(1 + H) - 4.*(1 + H).*C.*sin(2.*phi))).*a11 ...
                + (1 + 3.*H).*(sin(2.*theta).^2).*(B - C.*sin(2.*phi)).*a13 + 4.*(2 - 2.*H + (1 + 3.*H)...
                .*(sin(theta).^2)).*(A + (sin(theta).^2).*(B - C.*sin(2.*phi))).*a33));
        end
        
        function Qval = Qphi1(~, theta, phi, Beta1, I1, C, H, a11, a13, a33)
            % Beta_n is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent)
            % I_n is Intensity/(1-SBR^{-1}) (can be wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Qval = (C.^2).*(sin(theta).^4).*(sin(2.*phi).^2).*(((7 + 5.*H + (1 + 3.*H).*cos(2.*theta)).*a11...
                - 2.*(1 + 3.*H).*(cos(theta).^2).*a13 + (-5 + H + (1 + 3.*H).*cos(2.*theta)).*a33).^2)...
                .*(((1./Beta1) - 1).^2).*I1;
        end
        
        function Qval = Qphi1NB(~, theta, phi, I1, C, H, a11, a13, a33)
            % I_n is Intensity (can be wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Qval = (C.^2).*(sin(theta).^4).*(sin(2.*phi).^2).*(((7 + 5.*H + (1 + 3.*H).*cos(2.*theta)).*a11...
                - 2.*(1 + 3.*H).*(cos(theta).^2).*a13 + (-5 + H + (1 + 3.*H).*cos(2.*theta)).*a33).^2)...
                .*I1;
        end 
        
        function Mval = Mphi1(~, theta, phi, Beta1, A, B, C, H, NF, a11, a13, a33)
            % Beta is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent)
            % I is Intensity/(1-SBR^{-1}) (can be wavelength dependent),
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms            
            Mval = NF.*((4.*(2.*(A + B).*(3 + H) - (A + 3.*A.*H + 4.*B.*(1 + H) + ...
                4.*(1 + H).*C.*cos(2.*phi)).*(sin(theta).^2)...
                + (1 + 3.*H).*(B + C.*cos(2.*phi)).*(sin(theta).^4)).*a11 + (1 + 3.*H).*(B + C.*cos(2.*phi))...
                .*(sin(2.*theta).^2).*a13 + 4.*(2 - 2.*H + (1 + 3.*H).*(sin(theta).^2)).*(A + (B + C.*cos(2.*phi))...
                .*(sin(theta).^2)).*a33).*(1 - (1./Beta1)) + ((16.*NF)./Beta1));
        end
        
        function Mval = Mphi1NB(~, theta, phi, A, B, C, H, NF, a11, a13, a33)
            % I is Intensity (can be wavelength dependent),
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms            
            Mval = NF.*((4.*(2.*(A + B).*(3 + H) - (A + 3.*A.*H + 4.*B.*(1 + H) + ...
                4.*(1 + H).*C.*cos(2.*phi)).*(sin(theta).^2)...
                + (1 + 3.*H).*(B + C.*cos(2.*phi)).*(sin(theta).^4)).*a11 + (1 + 3.*H).*(B + C.*cos(2.*phi))...
                .*(sin(2.*theta).^2).*a13 + 4.*(2 - 2.*H + (1 + 3.*H).*(sin(theta).^2)).*(A + (B + C.*cos(2.*phi))...
                .*(sin(theta).^2)).*a33));
        end
        
        function Qval = Qphi2(~, theta, phi, Beta2, I2, C, H, a11, a13, a33)
            % Beta_n is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent)
            % I_n is Intensity/(1-SBR^{-1}) (can be wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Qval = (C.^2).*(sin(theta).^4).*(cos(2.*phi).^2).*(((7 + 5.*H + (1 + 3.*H).*cos(2.*theta)).*a11...
                - 2.*(1 + 3.*H).*(cos(theta).^2).*a13 + (-5 + H + (1 + 3.*H).*cos(2.*theta)).*a33).^2)...
                .*(((1./Beta2) - 1).^2).*I2;
        end
        
        function Qval = Qphi2NB(~, theta, phi, I2, C, H, a11, a13, a33)
            % I_n is Intensity (can be wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Qval = (C.^2).*(sin(theta).^4).*(cos(2.*phi).^2).*(((7 + 5.*H + (1 + 3.*H).*cos(2.*theta)).*a11...
                - 2.*(1 + 3.*H).*(cos(theta).^2).*a13 + (-5 + H + (1 + 3.*H).*cos(2.*theta)).*a33).^2)...
                .*I2;
        end

        
        function Mval = Mphi2(~, theta, phi, Beta2, A, B, C, H, NF, a11, a13, a33)
            % time is time taken to observe photons
            % Beta is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent)
            % I is Intensity/(1-SBR^{-1}) (can be wavelength dependent),
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Mval = NF.*((4.*(2.*(A + B).*(3 + H) + (1 + 3.*H).*(sin(theta).^4).*(B + C.*sin(2.*phi)) ...
                - (sin(theta).^2).*(A + 3.*A.*H + 4.*B.*(1 + H) + 4.*(1 + H).*C.*sin(2.*phi))).*a11 + ...
                (1 + 3.*H).*(sin(2.*theta).^2).*(B + C.*sin(2.*phi)).*a13 + 4.*(2 - 2.*H + (1 + 3.*H).*(sin(theta).^2))...
                .*(A + (sin(theta).^2).*(B + C.*sin(2.*phi))).*a33).*(1 - (1./Beta2)) + ((16.*NF)./Beta2));
        end
        
        function Mval = Mphi2NB(~, theta, phi, A, B, C, H, NF, a11, a13, a33)
            % I is Intensity (can be wavelength dependent),
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Mval = NF.*((4.*(2.*(A + B).*(3 + H) + (1 + 3.*H).*(sin(theta).^4).*(B + C.*sin(2.*phi)) ...
                - (sin(theta).^2).*(A + 3.*A.*H + 4.*B.*(1 + H) + 4.*(1 + H).*C.*sin(2.*phi))).*a11 + ...
                (1 + 3.*H).*(sin(2.*theta).^2).*(B + C.*sin(2.*phi)).*a13 + 4.*(2 - 2.*H + (1 + 3.*H).*(sin(theta).^2))...
                .*(A + (sin(theta).^2).*(B + C.*sin(2.*phi))).*a33));
        end
        
        function Qval = Qphi3(~, theta, phi, Beta3, I3, C, H, a11, a13, a33)
            % Beta_n is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent)
            % I_n is Intensity/(1-SBR^{-1}) (can be wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Qval = (C.^2).*(sin(theta).^4).*(sin(2.*phi).^2).*(((7 + 5.*H + (1 + 3.*H).*cos(2.*theta)).*a11...
                - 2.*(1 + 3.*H).*(cos(theta).^2).*a13 + (-5 + H + (1 + 3.*H).*cos(2.*theta)).*a33).^2)...
                .*(((1./Beta3) - 1).^2).*I3;
        end
        
        function Qval = Qphi3NB(~, theta, phi, I3, C, H, a11, a13, a33)
            % I_n is Intensity (can be wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Qval = (C.^2).*(sin(theta).^4).*(sin(2.*phi).^2).*(((7 + 5.*H + (1 + 3.*H).*cos(2.*theta)).*a11...
                - 2.*(1 + 3.*H).*(cos(theta).^2).*a13 + (-5 + H + (1 + 3.*H).*cos(2.*theta)).*a33).^2)...
                .*I3;
        end
        
        function Mval = Mphi3(~, theta, phi, Beta3, A, B, C, H, NF, a11, a13, a33)
            % Beta is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent)
            % I is Intensity/(1-SBR^{-1}) (can be wavelength dependent),
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Mval = NF.*((4.*(2.*(A + B).*(3 + H) - (A + 3.*A.*H + 4.*B.*(1 + H) - ...
                4.*(1 + H).*C.*cos(2.*phi)).*(sin(theta).^2)...
                + (1 + 3.*H).*(B - C.*cos(2.*phi)).*(sin(theta).^4)).*a11 + (1 + 3.*H).*(B - C.*cos(2.*phi))...
                .*(sin(2.*theta).^2).*a13 + 4.*(2 - 2.*H + (1 + 3.*H).*(sin(theta).^2)).*(A + (B - C.*cos(2.*phi))...
                .*(sin(theta).^2)).*a33).*(1 - (1./Beta3)) + ((16.*NF)./Beta3));
        end
        
        function Mval = Mphi3NB(~, theta, phi, A, B, C, H, NF, a11, a13, a33)
            % I is Intensity (can be wavelength dependent),
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Mval = NF.*((4.*(2.*(A + B).*(3 + H) - (A + 3.*A.*H + 4.*B.*(1 + H) - ...
                4.*(1 + H).*C.*cos(2.*phi)).*(sin(theta).^2)...
                + (1 + 3.*H).*(B - C.*cos(2.*phi)).*(sin(theta).^4)).*a11 + (1 + 3.*H).*(B - C.*cos(2.*phi))...
                .*(sin(2.*theta).^2).*a13 + 4.*(2 - 2.*H + (1 + 3.*H).*(sin(theta).^2)).*(A + (B - C.*cos(2.*phi))...
                .*(sin(theta).^2)).*a33));
        end
        
        function Qval = Qphi4(~, theta, phi, Beta4, I4, C, H, a11, a13, a33)
            % Beta_n is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent)
            % I_n is Intensity/(1-SBR^{-1}) (can be wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Qval = (C.^2).*(sin(theta).^4).*(cos(2.*phi).^2).*(((7 + 5.*H + (1 + 3.*H).*cos(2.*theta)).*a11...
                - 2.*(1 + 3.*H).*(cos(theta).^2).*a13 + (-5 + H + (1 + 3.*H).*cos(2.*theta)).*a33).^2)...
                .*(((1./Beta4) - 1).^2).*I4;
        end
        
        function Qval = Qphi4NB(~, theta, phi, I4, C, H, a11, a13, a33)
            % I_n is Intensity (can be wavelength dependent)
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Qval = (C.^2).*(sin(theta).^4).*(cos(2.*phi).^2).*(((7 + 5.*H + (1 + 3.*H).*cos(2.*theta)).*a11...
                - 2.*(1 + 3.*H).*(cos(theta).^2).*a13 + (-5 + H + (1 + 3.*H).*cos(2.*theta)).*a33).^2)...
                .*I4;
        end
        
        function Mval = Mphi4(~, theta, phi, Beta4, A, B, C, H, NF, a11, a13, a33)
            % time is time taken to observe photons
            % Beta is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent)
            % I is Intensity/(1-SBR^{-1}) (can be wavelength dependent),
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Mval = NF.*((4.*(2.*(A + B).*(3 + H) + (1 + 3.*H).*(sin(theta).^4).*(B - C.*sin(2.*phi)) ...
                - (sin(theta).^2).*(A + 3.*A.*H + 4.*B.*(1 + H) - 4.*(1 + H).*C.*sin(2.*phi))).*a11 + ...
                (1 + 3.*H).*(sin(2.*theta).^2).*(B - C.*sin(2.*phi)).*a13 + 4.*(2 - 2.*H + (1 + 3.*H).*(sin(theta).^2))...
                .*(A + (sin(theta).^2).*(B - C.*sin(2.*phi))).*a33).*(1 - (1./Beta4)) + ((16.*NF)./Beta4));
        end
        
        function Mval = Mphi4NB(~, theta, phi, A, B, C, H, NF, a11, a13, a33)
            % time is time taken to observe photons
            % Beta is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent)
            % I is Intensity/(1-SBR^{-1}) (can be wavelength dependent),
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            
            Mval = NF.*((4.*(2.*(A + B).*(3 + H) + (1 + 3.*H).*(sin(theta).^4).*(B - C.*sin(2.*phi)) ...
                - (sin(theta).^2).*(A + 3.*A.*H + 4.*B.*(1 + H) - 4.*(1 + H).*C.*sin(2.*phi))).*a11 + ...
                (1 + 3.*H).*(sin(2.*theta).^2).*(B - C.*sin(2.*phi)).*a13 + 4.*(2 - 2.*H + (1 + 3.*H).*(sin(theta).^2))...
                .*(A + (sin(theta).^2).*(B - C.*sin(2.*phi))).*a33));
        end
        
        function FisherInfoM = FisherMatrixScatter(obj, time, theta, phi, I, Beta,...
                A, B, C, H, NF, a11, a13, a33, n_detectors)
            % time is time taken to observe photons
            % Beta is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent), should be n_detector x wavelength
            % matrix
            % I is Intensity/(1-SBR^{-1}) (can be wavelength dependent),
            % should be n_detector x wavelength matrix
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms

            if size(I, 2) == 1 % if intensity is a simple scalar
                if length(A) ~= 1 % if we have multiple wavelengths
                    I = I./length(A); % spread intensity across multiple wavelengths
                end
            end
            
            if n_detectors == 3
                I1 = I(1, :); Beta1 = Beta(1, :);
                I2 = I(2, :); Beta2 = Beta(2, :);
                I3 = I(3, :); Beta3 = Beta(3, :);
                I4 = 0; Beta4 = 0;
            else
                I1 = I(1, :); Beta1 = Beta(1, :);
                I2 = I(2, :); Beta2 = Beta(2, :);
                I3 = I(3, :); Beta3 = Beta(3, :);
                I4 = I(4, :); Beta4 = Beta(4, :);
            end

            FisherInfoM = zeros(2,2);
            
            % get Q and M terms for thetatheta
            Qtheta1 = obj.Qthetatheta1(theta, phi, Beta1, I1, A, B, C, H, a11, a13, a33);
            Mtheta1 = obj.Mthetatheta1(theta, phi, Beta1, A, B, C, H, NF, a11, a13, a33);
            Qtheta2 = obj.Qthetatheta2(theta, phi, Beta2, I2, A, B, C, H, a11, a13, a33);
            Mtheta2 = obj.Mthetatheta2(theta, phi, Beta2, A, B, C, H, NF, a11, a13, a33);
            Qtheta3 = obj.Qthetatheta3(theta, phi, Beta3, I3, A, B, C, H, a11, a13, a33);
            Mtheta3 = obj.Mthetatheta3(theta, phi, Beta3, A, B, C, H, NF, a11, a13, a33);
            if n_detectors == 4
                Qtheta4 = obj.Qthetatheta4(theta, phi, Beta4, I4, A, B, C, H, a11, a13, a33);
                Mtheta4 = obj.Mthetatheta4(theta, phi, Beta4, A, B, C, H, NF, a11, a13, a33);
            end
            
            if n_detectors == 3
                FisherInfoM(1,1) = time.*(sum((Qtheta1./Mtheta1)) + sum((Qtheta2./Mtheta2)) +...
                    sum((Qtheta3./Mtheta3))); % calculate theta theta term
            elseif n_detectors == 4
                FisherInfoM(1,1) = time.*(sum((Qtheta1./Mtheta1)) + sum((Qtheta2./Mtheta2)) +...
                    sum((Qtheta3./Mtheta3)) + sum((Qtheta4./Mtheta4))); % calculate theta theta term
            end
            Qthetaphi1 = obj.Qthetaphi1(theta, phi, Beta1, I1, A, B, C, H, a11, a13, a33);
            Mthetaphi1 = obj.Mthetaphi1(theta, phi, Beta1, A, B, C, H, NF, a11, a13, a33);
            Qthetaphi2 = obj.Qthetaphi2(theta, phi, Beta2, I2, A, B, C, H, a11, a13, a33);
            Mthetaphi2 = obj.Mthetaphi2(theta, phi, Beta2, A, B, C, H, NF, a11, a13, a33);
            Qthetaphi3 = obj.Qthetaphi3(theta, phi, Beta3, I3, A, B, C, H, a11, a13, a33);
            Mthetaphi3 = obj.Mthetaphi3(theta, phi, Beta3, A, B, C, H, NF, a11, a13, a33);
            if n_detectors == 4
                Qthetaphi4 = obj.Qthetaphi4(theta, phi, Beta4, I4, A, B, C, H, a11, a13, a33);
                Mthetaphi4 = obj.Mthetaphi4(theta, phi, Beta4, A, B, C, H, NF, a11, a13, a33);
            end
            
            if n_detectors == 3
                FisherInfoM(1,2) = time.*(sum((Qthetaphi1./Mthetaphi1)) - sum((Qthetaphi2./Mthetaphi2)) -...
                    sum((Qthetaphi3./Mthetaphi3))); % calculate theta phi term
                FisherInfoM(2,1) = FisherInfoM(1,2);
            elseif n_detectors == 4
                FisherInfoM(1,2) = time.*(sum((Qthetaphi1./Mthetaphi1)) - sum((Qthetaphi2./Mthetaphi2)) -...
                    sum((Qthetaphi3./Mthetaphi3)) + sum((Qthetaphi4./Mthetaphi4))); % calculate theta phi term
                FisherInfoM(2,1) = FisherInfoM(1,2);
            end
            
            Qphi1 = obj.Qphi1(theta, phi, Beta1, I1, C, H, a11, a13, a33);
            Mphi1 = obj.Mphi1(theta, phi, Beta1, A, B, C, H, NF, a11, a13, a33);
            Qphi2 = obj.Qphi2(theta, phi, Beta2, I2, C, H, a11, a13, a33);
            Mphi2 = obj.Mphi2(theta, phi, Beta2, A, B, C, H, NF, a11, a13, a33);
            Qphi3 = obj.Qphi3(theta, phi, Beta3, I3, C, H, a11, a13, a33);
            Mphi3 = obj.Mphi3(theta, phi, Beta3, A, B, C, H, NF, a11, a13, a33);
            if n_detectors == 4
                Qphi4 = obj.Qphi4(theta, phi, Beta4, I4, C, H, a11, a13, a33);
                Mphi4 = obj.Mphi4(theta, phi, Beta4, A, B, C, H, NF, a11, a13, a33);
            end

            if n_detectors == 3
                FisherInfoM(2,2) = time.*(sum((Qphi1./Mphi1)) + sum((Qphi2./Mphi2)) +...
                    sum((Qphi3./Mphi3))); % calculate phi phi term
            elseif n_detectors == 4
                FisherInfoM(2,2) = time.*(sum((Qphi1./Mphi1)) + sum((Qphi2./Mphi2)) +...
                    sum((Qphi3./Mphi3)) + sum((Qphi4./Mphi4))); % calculate phi phi term
            end

        end
        
        function FisherInfoM = FisherMatrixScatterNB(obj, time, theta, phi, I,...
                A, B, C, H, NF, a11, a13, a33, n_detectors)
            % FIM for no background case
            % time is time taken to observe photons
            % I is Intensity (can be wavelength dependent),
            % should be n_detector x wavelength matrix
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms

            if size(I, 2) == 1 % if intensity is a simple scalar
                if length(A) ~= 1 % if we have multiple wavelengths
                    I = I./length(A); % spread intensity across multiple wavelengths
                end
            end
            
            if n_detectors == 3
                I1 = I(1, :);
                I2 = I(2, :);
                I3 = I(3, :);
                I4 = 0;
            else
                I1 = I(1, :);
                I2 = I(2, :);
                I3 = I(3, :);
                I4 = I(4, :);
            end

            FisherInfoM = zeros(2,2);
            
            % get Q and M terms for thetatheta
            Qtheta1 = obj.Qthetatheta1NB(theta, phi, I1, A, B, C, H, a11, a13, a33);
            Mtheta1 = obj.Mthetatheta1NB(theta, phi, A, B, C, H, NF, a11, a13, a33);
            Qtheta2 = obj.Qthetatheta2NB(theta, phi, I2, A, B, C, H, a11, a13, a33);
            Mtheta2 = obj.Mthetatheta2NB(theta, phi, A, B, C, H, NF, a11, a13, a33);
            Qtheta3 = obj.Qthetatheta3NB(theta, phi, I3, A, B, C, H, a11, a13, a33);
            Mtheta3 = obj.Mthetatheta3NB(theta, phi, A, B, C, H, NF, a11, a13, a33);
            if n_detectors == 4
                Qtheta4 = obj.Qthetatheta4NB(theta, phi, I4, A, B, C, H, a11, a13, a33);
                Mtheta4 = obj.Mthetatheta4NB(theta, phi, A, B, C, H, NF, a11, a13, a33);
            end
            
            if n_detectors == 3
                FisherInfoM(1,1) = time.*(sum((Qtheta1./Mtheta1)) + sum((Qtheta2./Mtheta2)) +...
                    sum((Qtheta3./Mtheta3))); % calculate theta theta term
            elseif n_detectors == 4
                FisherInfoM(1,1) = time.*(sum((Qtheta1./Mtheta1)) + sum((Qtheta2./Mtheta2)) +...
                    sum((Qtheta3./Mtheta3)) + sum((Qtheta4./Mtheta4))); % calculate theta theta term
            end
            Qthetaphi1 = obj.Qthetaphi1NB(theta, phi, I1, A, B, C, H, a11, a13, a33);
            Mthetaphi1 = obj.Mthetaphi1NB(theta, phi, A, B, C, H, NF, a11, a13, a33);
            Qthetaphi2 = obj.Qthetaphi2NB(theta, phi, I2, A, B, C, H, a11, a13, a33);
            Mthetaphi2 = obj.Mthetaphi2NB(theta, phi, A, B, C, H, NF, a11, a13, a33);
            Qthetaphi3 = obj.Qthetaphi3NB(theta, phi, I3, A, B, C, H, a11, a13, a33);
            Mthetaphi3 = obj.Mthetaphi3NB(theta, phi, A, B, C, H, NF, a11, a13, a33);
            if n_detectors == 4
                Qthetaphi4 = obj.Qthetaphi4NB(theta, phi, I4, A, B, C, H, a11, a13, a33);
                Mthetaphi4 = obj.Mthetaphi4NB(theta, phi, A, B, C, H, NF, a11, a13, a33);
            end
            
            if n_detectors == 3
                FisherInfoM(1,2) = time.*(sum((Qthetaphi1./Mthetaphi1)) - sum((Qthetaphi2./Mthetaphi2)) -...
                    sum((Qthetaphi3./Mthetaphi3))); % calculate theta phi term
                FisherInfoM(2,1) = FisherInfoM(1,2);
            elseif n_detectors == 4
                FisherInfoM(1,2) = time.*(sum((Qthetaphi1./Mthetaphi1)) - sum((Qthetaphi2./Mthetaphi2)) -...
                    sum((Qthetaphi3./Mthetaphi3)) + sum((Qthetaphi4./Mthetaphi4))); % calculate theta phi term
                FisherInfoM(2,1) = FisherInfoM(1,2);
            end
            
            Qphi1 = obj.Qphi1NB(theta, phi, I1, C, H, a11, a13, a33);
            Mphi1 = obj.Mphi1NB(theta, phi, A, B, C, H, NF, a11, a13, a33);
            Qphi2 = obj.Qphi2NB(theta, phi, I2, C, H, a11, a13, a33);
            Mphi2 = obj.Mphi2NB(theta, phi, A, B, C, H, NF, a11, a13, a33);
            Qphi3 = obj.Qphi3NB(theta, phi, I3, C, H, a11, a13, a33);
            Mphi3 = obj.Mphi3NB(theta, phi, A, B, C, H, NF, a11, a13, a33);
            if n_detectors == 4
                Qphi4 = obj.Qphi4NB(theta, phi, I4, C, H, a11, a13, a33);
                Mphi4 = obj.Mphi4NB(theta, phi, A, B, C, H, NF, a11, a13, a33);
            end

            if n_detectors == 3
                FisherInfoM(2,2) = time.*(sum((Qphi1./Mphi1)) + sum((Qphi2./Mphi2)) +...
                    sum((Qphi3./Mphi3))); % calculate phi phi term
            elseif n_detectors == 4
                FisherInfoM(2,2) = time.*(sum((Qphi1./Mphi1)) + sum((Qphi2./Mphi2)) +...
                    sum((Qphi3./Mphi3)) + sum((Qphi4./Mphi4))); % calculate phi phi term
            end

        end
        
        function FisherInfoM = FisherMatrixFluo(obj, time, theta, phi, I, Beta,...
                A, B, C, n_detectors)
            % time is time taken to observe photons
            % Beta is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent), should be n_detector x wavelength
            % matrix
            % I is Intensity/(1-SBR^{-1}) (can be wavelength dependent),
            % should be n_detector x wavelength matrix
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            FisherInfoM = zeros(2,2);
            FisherInfoM(1,1) = obj.Jthetathetaf(theta, phi, time, Beta, I,...
                A, B, C, n_detectors);
            FisherInfoM(1,2) = obj.Jthetaphif(theta, phi, time, Beta, I,...
                A, B, C, n_detectors);
            FisherInfoM(2,1) = FisherInfoM(1,2);
            FisherInfoM(2,2) = obj.Jphiphif(theta, phi, time, Beta, I,...
                A, B, C, n_detectors);
        end
        
        function FisherInfoM = FisherMatrixFluoNB(obj, time, theta, phi, I,...
                A, B, C, n_detectors)
            % FIM for no background case
            % time is time taken to observe photons
            % I is Intensity (can be wavelength dependent),
            % should be n_detector x wavelength matrix
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            FisherInfoM = zeros(2,2);
            FisherInfoM(1,1) = obj.JthetathetafNB(theta, phi, time, I,...
                A, B, C, n_detectors);
            FisherInfoM(1,2) = obj.JthetaphifNB(theta, phi, time, I,...
                A, B, C, n_detectors);
            FisherInfoM(2,1) = FisherInfoM(1,2);
            FisherInfoM(2,2) = obj.JphiphifNB(theta, phi, time, I,...
                A, B, C, n_detectors);
        end
        
        function [CramerRaoMTheta, CramerRaoMPhi] = CramerRaoScatter(obj, time, thetas, phis, I, Beta,...
                wavelengths, T, cv, NAObj, NACond, a11, a13, a33, n_detectors, thetadepf)
            % time is time taken to observe photons
            % Beta is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent), should be n_detector x wavelength
            % matrix
            % I is Intensity/(1-SBR^{-1}) (can be wavelength dependent),
            % should be n_detector x wavelength matrix
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            % thetadepf is an optional argument that weights the
            % intensity such that over a freely-diffusing object the
            % intensity would average to I
            CramerRaoMTheta = zeros(length(thetas), length(phis));
            CramerRaoMPhi = zeros(length(thetas), length(phis));
            n0 = obj.n_m(wavelengths, T, cv);
            [~, ~, A, B, C, H] = obj.InstrResp(NACond, NAObj, n0);
            for i = 1:length(thetas)
                NF = obj.NormFactor(thetas(i), a11, a13, a33, A, B, H);
                if nargin>=16
                    Ival = I.*thetadepf(thetas(i));
                else
                    Ival = I;
                end
                for j = 1:length(phis)
                    J = obj.FisherMatrixScatter(time, thetas(i), phis(j), Ival, Beta,...
                A, B, C, H, NF, a11, a13, a33, n_detectors);
                    InvJ = inv(J);
                    CramerRaoMTheta(i,j) = real(sqrt(InvJ(1,1)));
                    CramerRaoMPhi(i,j) = real(sqrt(InvJ(2,2))); 
                end
            end
            % cannot be "worse" than pi/2 and pi respectively
            CramerRaoMTheta(CramerRaoMTheta > pi./2) = pi./2;
            CramerRaoMPhi(CramerRaoMPhi > pi) = pi;
            CramerRaoMTheta(isnan(CramerRaoMTheta)) = pi./2;
            CramerRaoMPhi(isnan(CramerRaoMPhi)) = pi;
        end
        
        function [CramerRaoMTheta, CramerRaoMPhi] = CramerRaoScatterOOI(obj, time, thetas, phis, I, Beta,...
                wavelengths, T, cv, NAObj, NACond, a11, a13, a33, n_detectors, thetadepf)
            % time is time taken to observe photons
            % Beta is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent), should be n_detector x wavelength
            % matrix
            % I is Intensity/(1-SBR^{-1}) (can be wavelength dependent),
            % should be n_detector x wavelength matrix
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            % thetadepf is an optional argument that weights the
            % intensity such that over a freely-diffusing object the
            % intensity would average to I
            CramerRaoMTheta = zeros(length(thetas), length(phis));
            CramerRaoMPhi = zeros(length(thetas), length(phis));
            n0 = obj.n_m(wavelengths, T, cv);
            [~, ~, A, B, C, H] = obj.InstrResp(NACond, NAObj, n0);
            for i = 1:length(thetas)
                NF = obj.NormFactor(thetas(i), a11, a13, a33, A, B, H);
                if nargin>=16
                    Ival = I.*thetadepf(thetas(i));
                else
                    Ival = I;
                end
                for j = 1:length(phis)
                    J = obj.FisherMatrixScatter(time, thetas(i), phis(j), Ival, Beta,...
                A, B, C, H, NF, a11, a13, a33, n_detectors);
                    InvJ = inv(J);
                    CramerRaoMTheta(i,j) = real(sqrt(InvJ(1,1)));
                    CramerRaoMPhi(i,j) = real(sqrt(InvJ(2,2))); 
                end
            end
        end

        function CramerRaoMTheta = CramerRaoScatterThetErr(obj, time, thetas, phis, I, Beta,...
                wavelengths, T, cv, NAObj, NACond, a11, a13, a33, n_detectors, thetadepf)
            % time is time taken to observe photons
            % Beta is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent), should be n_detector x wavelength
            % matrix
            % I is Intensity/(1-SBR^{-1}) (can be wavelength dependent),
            % should be n_detector x wavelength matrix
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            % thetadepf is an optional argument that weights the
            % intensity such that over a freely-diffusing object the
            % intensity would average to I
            CramerRaoMTheta = zeros(1, length(thetas));
            n0 = obj.n_m(wavelengths, T, cv);
            [~, ~, A, B, C, H] = obj.InstrResp(NACond, NAObj, n0);
            for i = 1:length(thetas)
                NF = obj.NormFactor(thetas(i), a11, a13, a33, A, B, H);
                if nargin>=16
                    Ival = I.*thetadepf(thetas(i));
                else
                    Ival = I;
                end
                J = obj.FisherMatrixScatter(time, thetas(i), phis(i), Ival, Beta,...
            A, B, C, H, NF, a11, a13, a33, n_detectors);
                InvJ = inv(J);
                CramerRaoMTheta(1,i) = real(sqrt(InvJ(1,1)));
            end
            % cannot be "worse" than pi/2 and pi respectively
            CramerRaoMTheta(CramerRaoMTheta > pi./2) = pi./2;
            CramerRaoMTheta(isnan(CramerRaoMTheta)) = pi./2;
        end
        
        function CramerRaoMPhi = CramerRaoScatterPhiErr(obj, time, thetas, phis, I, Beta,...
                wavelengths, T, cv, NAObj, NACond, a11, a13, a33, n_detectors, thetadepf)
            % time is time taken to observe photons
            % Beta is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent), should be n_detector x wavelength
            % matrix
            % I is Intensity/(1-SBR^{-1}) (can be wavelength dependent),
            % should be n_detector x wavelength matrix
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            % thetadepf is an optional argument that weights the
            % intensity such that over a freely-diffusing object the
            % intensity would average to I
            CramerRaoMPhi = zeros(1, length(phis));
            n0 = obj.n_m(wavelengths, T, cv);
            [~, ~, A, B, C, H] = obj.InstrResp(NACond, NAObj, n0);
            for i = 1:length(thetas)
                NF = obj.NormFactor(thetas(i), a11, a13, a33, A, B, H);
                if nargin>=16
                    Ival = I.*thetadepf(thetas(i));
                else
                    Ival = I;
                end
                J = obj.FisherMatrixScatter(time, thetas(i), phis(i), Ival, Beta,...
            A, B, C, H, NF, a11, a13, a33, n_detectors);
                InvJ = inv(J);
                CramerRaoMPhi(1,i) = real(sqrt(InvJ(2,2))); 
            end
            % cannot be "worse" than pi/2 and pi respectively
            CramerRaoMPhi(CramerRaoMPhi > pi) = pi;
            CramerRaoMPhi(isnan(CramerRaoMPhi)) = pi;
        end

        function [CramerRaoMTheta, CramerRaoMPhi] = CramerRaoScatterNB(obj, time, thetas, phis, I,...
                wavelengths, T, cv, NAObj, NACond, a11, a13, a33, n_detectors, thetadepf)
            % CRLB for no-background case
            % time is time taken to observe photons
            % I is Intensity (can be wavelength dependent),
            % should be n_detector x wavelength matrix
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            % thetadepf is an optional argument that weights the
            % intensity such that over a freely-diffusing object the
            % intensity would average to I            
            CramerRaoMTheta = zeros(length(thetas), length(phis));
            CramerRaoMPhi = zeros(length(thetas), length(phis));
            n0 = obj.n_m(wavelengths, T, cv);
            [~, ~, A, B, C, H] = obj.InstrResp(NACond, NAObj, n0);
            for i = 1:length(thetas)
                NF = obj.NormFactor(thetas(i), a11, a13, a33, A, B, H);
                if nargin>=15
                    Ival = I.*thetadepf(thetas(i));
                else
                    Ival = I;
                end
                for j = 1:length(phis)
                    J = obj.FisherMatrixScatterNB(time, thetas(i), phis(j), Ival,...
                A, B, C, H, NF, a11, a13, a33, n_detectors);
                    InvJ = inv(J);
                    CramerRaoMTheta(i,j) = real(sqrt(InvJ(1,1)));
                    CramerRaoMPhi(i,j) = real(sqrt(InvJ(2,2))); 
                end
            end
            % cannot be "worse" than pi/2 and pi respectively
            CramerRaoMTheta(CramerRaoMTheta > pi./2) = pi./2;
            CramerRaoMPhi(CramerRaoMPhi > pi) = pi;
            CramerRaoMTheta(isnan(CramerRaoMTheta)) = pi./2;
            CramerRaoMPhi(isnan(CramerRaoMPhi)) = pi;
        end
        
        function [CramerRaoMTheta, CramerRaoMPhi] = CramerRaoScatterNBOOI(obj, time, thetas, phis, I,...
                wavelengths, T, cv, NAObj, NACond, a11, a13, a33, n_detectors, thetadepf)
            % CRLB for no-background case
            % time is time taken to observe photons
            % I is Intensity (can be wavelength dependent),
            % should be n_detector x wavelength matrix
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % a11, a13 and a33 are polarizability tensor * lamp intensity
            % terms
            % thetadepf is an optional argument that weights the
            % intensity such that over a freely-diffusing object the
            % intensity would average to I            
            CramerRaoMTheta = zeros(length(thetas), length(phis));
            CramerRaoMPhi = zeros(length(thetas), length(phis));
            n0 = obj.n_m(wavelengths, T, cv);
            [~, ~, A, B, C, H] = obj.InstrResp(NACond, NAObj, n0);
            for i = 1:length(thetas)
                NF = obj.NormFactor(thetas(i), a11, a13, a33, A, B, H);
                if nargin>=15
                    Ival = I.*thetadepf(thetas(i));
                else
                    Ival = I;
                end
                for j = 1:length(phis)
                    J = obj.FisherMatrixScatterNB(time, thetas(i), phis(j), Ival,...
                A, B, C, H, NF, a11, a13, a33, n_detectors);
                    InvJ = inv(J);
                    CramerRaoMTheta(i,j) = real(sqrt(InvJ(1,1)));
                    CramerRaoMPhi(i,j) = real(sqrt(InvJ(2,2))); 
                end
            end
        end

        
        function [CramerRaoMTheta, CramerRaoMPhi] = CramerRaoFluoNB(obj, time, thetas, phis, I,...
                 wavelengths, T, cv, NAObj, n_detectors, thetadepf)
            % CRLB for non-background case
            % time is time taken to observe photons
            % I is Intensity (can be wavelength dependent),
            % should be n_detector x wavelength matrix
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % thetadepf is an optional argument that weights the
            % intensity such that over a freely-diffusing object the
            % intensity would average to I (computed separately)
            CramerRaoMTheta = zeros(length(thetas), length(phis));
            CramerRaoMPhi = zeros(length(thetas), length(phis));
            n0 = obj.n_m(wavelengths, T, cv);
            [~, ~, A, B, C, ~] = obj.InstrResp(1, NAObj, n0);
            for i = 1:length(thetas)
                if nargin>=11
                    Ival = I.*thetadepf(thetas(i));
                else
                    Ival = I;
                end
                for j = 1:length(phis)
                    J = obj.FisherMatrixFluoNB(time, thetas(i), phis(j), Ival,...
                A, B, C, n_detectors);
                    InvJ = inv(J);
                    CramerRaoMTheta(i,j) = sqrt(InvJ(1,1));
                    CramerRaoMPhi(i,j) = sqrt(InvJ(2,2)); 
                end
            end
            % cannot be "worse" than pi/2 and pi respectively
            CramerRaoMTheta(CramerRaoMTheta > pi./2) = pi./2;
            CramerRaoMPhi(CramerRaoMPhi > pi) = pi;
            CramerRaoMTheta(isnan(CramerRaoMTheta)) = pi./2;
            CramerRaoMPhi(isnan(CramerRaoMPhi)) = pi;
        end
        
        function [CramerRaoMTheta, CramerRaoMPhi] = CramerRaoFluoNBOTI(obj, time, thetas, phis, I,...
                 wavelengths, T, cv, NAObj, n_detectors, thetadepf)
            % CRLB for non-background case
            % time is time taken to observe photons
            % I is Intensity (can be wavelength dependent),
            % should be n_detector x wavelength matrix
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % thetadepf is an optional argument that weights the
            % intensity such that over a freely-diffusing object the
            % intensity would average to I (computed separately)
            CramerRaoMTheta = zeros(length(thetas), length(phis));
            CramerRaoMPhi = zeros(length(thetas), length(phis));
            n0 = obj.n_m(wavelengths, T, cv);
            [~, ~, A, B, C, ~] = obj.InstrResp(1, NAObj, n0);
            for i = 1:length(thetas)
                if nargin>=11
                    Ival = I.*thetadepf(thetas(i));
                else
                    Ival = I;
                end
                for j = 1:length(phis)
                    J = obj.FisherMatrixFluoNB(time, thetas(i), phis(j), Ival,...
                A, B, C, n_detectors);
                    InvJ = inv(J);
                    CramerRaoMTheta(i,j) = sqrt(InvJ(1,1));
                    CramerRaoMPhi(i,j) = sqrt(InvJ(2,2)); 
                end
            end
        end

        
        function [CramerRaoMTheta, CramerRaoMPhi] = CramerRaoFluo(obj, time, thetas, phis, I, Beta,...
                 wavelengths, T, cv, NAObj, n_detectors, thetadepf)
            % time is time taken to observe photons
            % Beta is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent), should be n_detector x wavelength
            % matrix
            % I is Intensity/(1-SBR^{-1}) (can be wavelength dependent),
            % should be n_detector x wavelength matrix
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % thetadepf is an optional argument that weights the
            % intensity such that over a freely-diffusing object the
            % intensity would average to I (computed separately)
            CramerRaoMTheta = zeros(length(thetas), length(phis));
            CramerRaoMPhi = zeros(length(thetas), length(phis));
            n0 = obj.n_m(wavelengths, T, cv);
            [~, ~, A, B, C, ~] = obj.InstrResp(1, NAObj, n0);
            for i = 1:length(thetas)
                if nargin>=12
                    Ival = I.*thetadepf(thetas(i));
                else
                    Ival = I;
                end
                for j = 1:length(phis)
                    J = obj.FisherMatrixFluo(time, thetas(i), phis(j), Ival, Beta,...
                A, B, C, n_detectors);
                    InvJ = inv(J);
                    CramerRaoMTheta(i,j) = sqrt(InvJ(1,1));
                    CramerRaoMPhi(i,j) = sqrt(InvJ(2,2)); 
                end
            end
            % cannot be "worse" than pi/2 and pi respectively
            CramerRaoMTheta(CramerRaoMTheta > pi./2) = pi./2;
            CramerRaoMPhi(CramerRaoMPhi > pi) = pi;
            CramerRaoMTheta(isnan(CramerRaoMTheta)) = pi./2;
            CramerRaoMPhi(isnan(CramerRaoMPhi)) = pi;
        end

        function [CramerRaoMTheta, CramerRaoMPhi] = CramerRaoFluoOTI(obj, time, thetas, phis, I, Beta,...
                 wavelengths, T, cv, NAObj, n_detectors, thetadepf)
            % time is time taken to observe photons
            % Beta is SBR = (Intensity + Background )/ Background (can be
            % wavelength dependent), should be n_detector x wavelength
            % matrix
            % I is Intensity/(1-SBR^{-1}) (can be wavelength dependent),
            % should be n_detector x wavelength matrix
            % wavelengths are wavelengths that hit the detector
            % n_detectors is how many detectors used
            % if doing multiple wavelengths, give wavelengths as row vector
            % and multiple angles as column vector
            % thetadepf is an optional argument that weights the
            % intensity such that over a freely-diffusing object the
            % intensity would average to I (computed separately)
            CramerRaoMTheta = zeros(length(thetas), length(phis));
            CramerRaoMPhi = zeros(length(thetas), length(phis));
            n0 = obj.n_m(wavelengths, T, cv);
            [~, ~, A, B, C, ~] = obj.InstrResp(1, NAObj, n0);
            for i = 1:length(thetas)
                if nargin>=12
                    Ival = I.*thetadepf(thetas(i));
                else
                    Ival = I;
                end
                for j = 1:length(phis)
                    J = obj.FisherMatrixFluo(time, thetas(i), phis(j), Ival, Beta,...
                A, B, C, n_detectors);
                    InvJ = inv(J);
                    CramerRaoMTheta(i,j) = sqrt(InvJ(1,1));
                    CramerRaoMPhi(i,j) = sqrt(InvJ(2,2)); 
                end
            end
        end
        
    end
end

