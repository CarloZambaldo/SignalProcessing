function n = noise(time,amplitude)
        
    if nargin < 2
        amplitude = 1;
    end
    n = amplitude*rand(1, length(time));

end