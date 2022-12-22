function y = fitness_TopoBARTMAP(x, data, classes)

    nIndividuals = size(x, 1);
    y = zeros(nIndividuals, 1);

    parfor ix =1:nIndividuals
        y(ix) = run_TopoBARTMAP(x(ix, :), data, classes);
    end

end
