function[impulse] = impulse_response(gamma1,gamma3,periods)
%only considers a system of the form:
%X_t = gamma_1 X_{t-1} + gamma_3 + \eta_t , \eta_t ~NID(0,sigma^2).

numShocks=size(gamma3,2);
numVars=size(gamma3,1);

impulse=zeros(periods,numVars,numShocks);

    for ii=1:numShocks
    shock=zeros(numShocks,1);
    shock(ii)=1;



    for jj=1:periods
    impulse(jj,:,ii)=gamma1^(jj-1)*gamma3*shock;    
    end
   
    end


end