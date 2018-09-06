figure('Name','Simulated Variables');
for jj=1:numVar
    subplot(5,4,jj);
    plot(XX(jj,:));
end

figure('Name','learning-intercept coefficients');
for jj=1:numEndo
    subplot(5,3,jj)
    plot(learning_matrix(:,1,jj));
end

figure('Name','learning-lagged coefficients');
index=0;
for jj=1:backward_indices
    for ii=1:numBackward
        index=index+1;
       subplot(numBackward,numBackward,index);
       plot(learning_matrix(:,1+ii,jj));
    end
end

figure('Name','learning-shock coefficients');
index=0;
for jj=1:numEndo
    for ii=1:numShocks
       index=index+1;
       subplot(numEndo,numShocks,index);
       plot(learning_matrix(:,numBackward+1+ii,jj));
    end
end

        
figure('Name','largest eigenvalue');
subplot(2,1,1);
plot(largestEig);
subplot(2,1,2);
plot(pr_flag);

