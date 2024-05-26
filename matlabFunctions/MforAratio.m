function [roots] = MforAratio(Aratio,k)
    % Call this script as [M] = MforAratio(A/A_*,gamma);
    % first argument is A/A_*
    % second argument is specific heat ratio
    
    % M will be an array of dimension 2
    % M(1) will be the subsonic root
    % M(2) will be the supersonic root

    if (Aratio < 1)
        disp(['Area ratio is less than 1: ' num2str(Aratio)])
        roots(1) = 1;
        roots(2) = 1;
        return
    end

    % Find subsonic root
    m = 0.5; % First guess
    % Newton Rhapson
    for iter=1:100
        equation = Aratio*m - ((1+(k-1)/2*m^2)/((k+1)/2))^((k+1)/(k-1)/2);  % Equation 9.44 in Cengal (without typo)
        diffeqm = Aratio - (m*(((1/2*k - 1/2)*m^2 + 1)/(1/2*k + 1/2))^((1/2*(k + 1))/(k - 1) - 1)*(1/2*k - 1/2)*(k + 1))/((1/2*k + 1/2)*(k - 1));
        m = m - equation/diffeqm;
        if (abs(equation) < 1.0e-8)
            break
        end

        if (m < 0 || m > 1)
            % Bad things are happening.  
            % Let's start again in between 0 and 1.
            m = rand(1);
        end
    end
    if (iter > 99)
        disp(['trouble converging for subsonic root with Aratio ' num2str(Aratio)])
        roots(1) = nan;
    else
        roots(1) = m;
    end

    % Find supersonic root
    m = 2.0; % First guess
    % Newton Rhapson
    for iter=1:100
        equation = Aratio*m - ((1+(k-1)/2*m^2)/((k+1)/2))^((k+1)/(k-1)/2);  % Equation 9.44 in Cengal (without typo)
        diffeqm = Aratio - (m*(((1/2*k - 1/2)*m^2 + 1)/(1/2*k + 1/2))^((1/2*(k + 1))/(k - 1) - 1)*(1/2*k - 1/2)*(k + 1))/((1/2*k + 1/2)*(k - 1));

        m = m - equation/diffeqm;
        if (abs(equation) < 1.0e-8)
            break
        end
        if (m < 1 | m > 1000)
            % Bad things are happening.  
            % Let's start again in between 1 and 1000.
            m = rand(1)*999 +1;
        end
    end
    if (iter > 99)
        disp(['trouble converging for supersonic root with Aratio ' num2str(Aratio)])
        roots(2) = nan;
    else
        roots(2) = m;
    end      
end
