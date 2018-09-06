figure('Name','learning-intercept coefficients');
for jj=1:numEndo
    subplot(5,3,jj);
    plot(matrix_learning(:,jj,1));
end

figure('Name','learning-autocorrelation coefficients');
for jj=1:numEndo
    subplot(5,3,jj)
    plot(matrix_learning(:,jj,2));
end

figure('Name','learning-simulated endogenous variables');
for jj=1:numEndo
    subplot(5,3,jj)
    plot(XX(jj,:));
end

