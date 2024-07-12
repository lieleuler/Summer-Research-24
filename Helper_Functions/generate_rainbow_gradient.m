
function gradient = generate_rainbow_gradient(steps)
    if steps < 4
        if steps == 1
            gradient = [1, 0, 0];
        elseif steps == 2
            gradient = [1, 0, 0;
                        1, 0, 1];
        elseif steps == 3
            gradient = [1, 0, 0;
                        0, 1, 0;
                        0, 0, 1];
        elseif steps == 4
            gradient = [1, 0, 0;
                        0, 1, 0;
                        0, 0, 1;
                        0, 1, 1];
        end
        return
    end

    N = floor(steps/5);
    gradient = zeros(5*N, 3);

    % Red to Yellow
    for i = 0:(N-1)
        gradient(i + 1, 1) = 1;
        gradient(i + 1, 2) = i/N; % 0 to just before 1
        gradient(i + 1, 3) = 0;
    end
    % Yellow to Green
    for i = N:(2*N - 1)
        gradient(i + 1, 1) = (2*N - i)/N; % 1 to just before 0
        gradient(i + 1, 2) = 1;
        gradient(i + 1, 3) = 0;
    end
    % Yellow to Green
    for i = 2*N:(3*N - 1)
        gradient(i + 1, 1) = 0;
        gradient(i + 1, 2) = 1;
        gradient(i + 1, 3) = (i - 2*N)/N; % 0 to jump before 1
    end
    % Green to Blue
    for i = 3*N:(4*N - 1)
        gradient(i + 1, 1) = 0; 
        gradient(i + 1, 2) = (4*N - i)/N; % 1 to just before 0
        gradient(i + 1, 3) = 1;
    end
    % Blue to Purple
    for i = 4*N:(5*N - 1)
        gradient(i + 1, 1) = (i - 4*N)/N; % 0 to just before 1
        gradient(i + 1, 2) = 0; % 1 to just before
        gradient(i + 1, 3) = 1;
    end
    % Append purple to end
    for j = 1:(steps - 5*N)
        gradient(5*N + j, 1) = 1; 
        gradient(5*N + j, 2) = 0;
        gradient(5*N + j, 3) = 1;
    end
end

