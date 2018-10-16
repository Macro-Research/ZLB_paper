function parasim = prior_draws(nsimul,bounds)
a=0.4;b=0.25;
i=1;
parasim(:,1)=normrnd(a,b,[nsimul 1]);
outbound = ~((bounds(i,1) < parasim(:,i)).*(parasim(:,i) < bounds(i,2)));
outboundset = find(outbound);

            for j=1:1:sum(outbound)
                notvalid = 1;
                while notvalid
                    temp_para = normrnd(a,b,1,1);
                    
                    if ((bounds(i,1) < temp_para)&&(temp_para < bounds(i,2)))
                        parasim(outboundset(j),1) = temp_para;
                        notvalid = 0;
                    end
                end
            end

i=2;
mu=0.62;sigma2=0.25^2;b  = sigma2/mu;a  = mu/b;
parasim(:,2)=gamrnd(a,b,[nsimul 1]);
outbound = ~((bounds(i,1) < parasim(:,i)).*(parasim(:,i) < bounds(i,2)));
outboundset = find(outbound);

            for j=1:1:sum(outbound)
                notvalid = 1;
                while notvalid
                    temp_para = gamrnd(a,b,1,1);
                    
                    if ((bounds(i,1) < temp_para)&&(temp_para < bounds(i,2)))
                        parasim(outboundset(j),2) = temp_para;
                        notvalid = 0;
                    end
                end
            end



i=3;
mu=0.5;sigma2=0.25^2;b  = sigma2/mu;a  = mu/b;
parasim(:,3)=gamrnd(a,b,[nsimul 1]);
outbound = ~((bounds(i,1) < parasim(:,i)).*(parasim(:,i) < bounds(i,2)));
outboundset = find(outbound);

            for j=1:1:sum(outbound)
                notvalid = 1;
                while notvalid
                    temp_para = gamrnd(a,b,1,1);
                    
                    if ((bounds(i,1) < temp_para)&&(temp_para < bounds(i,2)))
                        parasim(outboundset(j),2) = temp_para;
                        notvalid = 0;
                    end
                end
            end



i=4;
mu=0.2;sigma2=0.1^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
parasim(:,4)= betarnd(a,b,[nsimul 1]);
outbound = ~((bounds(i,1) < parasim(:,i)).*(parasim(:,i) < bounds(i,2)));
outboundset = find(outbound);

            for j=1:1:sum(outbound)
                notvalid = 1;
                while notvalid
                    temp_para = betarnd(a,b,1,1);
                    
                    if ((bounds(i,1) < temp_para)&&(temp_para < bounds(i,2)))
                        parasim(outboundset(j),2) = temp_para;
                        notvalid = 0;
                    end
                end
            end



i=5;
mu=2;sigma2=0.5^2;b  = sigma2/mu;a  = mu/b;
parasim(:,5)=gamrnd(a,b,[nsimul 1]);
outbound = ~((bounds(i,1) < parasim(:,i)).*(parasim(:,i) < bounds(i,2)));
outboundset = find(outbound);

            for j=1:1:sum(outbound)
                notvalid = 1;
                while notvalid
                    temp_para = gamrnd(a,b,1,1);
                    
                    if ((bounds(i,1) < temp_para)&&(temp_para < bounds(i,2)))
                        parasim(outboundset(j),2) = temp_para;
                        notvalid = 0;
                    end
                end
            end

i=6;
mu=1.5;sigma2=0.25^2;b  = sigma2/mu;a  = mu/b;
parasim(:,6)=gamrnd(a,b,[nsimul 1]);
outbound = ~((bounds(i,1) < parasim(:,i)).*(parasim(:,i) < bounds(i,2)));
outboundset = find(outbound);

            for j=1:1:sum(outbound)
                notvalid = 1;
                while notvalid
                    temp_para = gamrnd(a,b,1,1);
                    
                    if ((bounds(i,1) < temp_para)&&(temp_para < bounds(i,2)))
                        parasim(outboundset(j),2) = temp_para;
                        notvalid = 0;
                    end
                end
            end

i=7;
mu=0.5;sigma2=0.25^2;b  = sigma2/mu;a  = mu/b;
parasim(:,7)=gamrnd(a,b,[nsimul 1]);
outbound = ~((bounds(i,1) < parasim(:,i)).*(parasim(:,i) < bounds(i,2)));
outboundset = find(outbound);

            for j=1:1:sum(outbound)
                notvalid = 1;
                while notvalid
                    temp_para = gamrnd(a,b,1,1);
                    
                    if ((bounds(i,1) < temp_para)&&(temp_para < bounds(i,2)))
                        parasim(outboundset(j),2) = temp_para;
                        notvalid = 0;
                    end
                end
            end

i=8;
mu=0.5;sigma2=0.2^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
parasim(:,8)= betarnd(a,b,[nsimul 1]);
outbound = ~((bounds(i,1) < parasim(:,i)).*(parasim(:,i) < bounds(i,2)));
outboundset = find(outbound);

            for j=1:1:sum(outbound)
                notvalid = 1;
                while notvalid
                    temp_para = betarnd(a,b,1,1);
                    
                    if ((bounds(i,1) < temp_para)&&(temp_para < bounds(i,2)))
                        parasim(outboundset(j),2) = temp_para;
                        notvalid = 0;
                    end
                end
            end

i=9;
mu=0.5;sigma2=0.2^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
parasim(:,9)= betarnd(a,b,[nsimul 1]);
outbound = ~((bounds(i,1) < parasim(:,i)).*(parasim(:,i) < bounds(i,2)));
outboundset = find(outbound);

            for j=1:1:sum(outbound)
                notvalid = 1;
                while notvalid
                    temp_para = betarnd(a,b,1,1);
                    
                    if ((bounds(i,1) < temp_para)&&(temp_para < bounds(i,2)))
                        parasim(outboundset(j),2) = temp_para;
                        notvalid = 0;
                    end
                end
            end

i=10;
parasim(:,10)= betarnd(a,b,[nsimul 1]);
outbound = ~((bounds(i,1) < parasim(:,i)).*(parasim(:,i) < bounds(i,2)));
outboundset = find(outbound);

            for j=1:1:sum(outbound)
                notvalid = 1;
                while notvalid
                    temp_para = betarnd(a,b,1,1);
                    
                    if ((bounds(i,1) < temp_para)&&(temp_para < bounds(i,2)))
                        parasim(outboundset(j),2) = temp_para;
                        notvalid = 0;
                    end
                end
            end

i=11;
a = 0.7;b=2;parasim(:,11) = sqrt( b*a^2./sum( (randn(b,nsimul)).^2 )' );
outbound = ~((bounds(i,1) < parasim(:,i)).*(parasim(:,i) < bounds(i,2)));
outboundset = find(outbound);

            for j=1:1:sum(outbound)
                notvalid = 1;
                while notvalid
                    temp_para = sqrt( b*a^2./sum( (randn(b,1)).^2 )' );
                    
                    if ((bounds(i,1) < temp_para)&&(temp_para < bounds(i,2)))
                        parasim(outboundset(j),2) = temp_para;
                        notvalid = 0;
                    end
                end
            end

i=12;
a = 0.3;b=2;parasim(:,12) = sqrt( b*a^2./sum( (randn(b,nsimul)).^2 )' );
outbound = ~((bounds(i,1) < parasim(:,i)).*(parasim(:,i) < bounds(i,2)));
outboundset = find(outbound);

            for j=1:1:sum(outbound)
                notvalid = 1;
                while notvalid
                    temp_para = sqrt( b*a^2./sum( (randn(b,1)).^2 )' );
                    
                    if ((bounds(i,1) < temp_para)&&(temp_para < bounds(i,2)))
                        parasim(outboundset(j),2) = temp_para;
                        notvalid = 0;
                    end
                end
            end

i=13;
a = 0.3;b=2;parasim(:,13) = sqrt( b*a^2./sum( (randn(b,nsimul)).^2 )' );
outbound = ~((bounds(i,1) < parasim(:,i)).*(parasim(:,i) < bounds(i,2)));
outboundset = find(outbound);

            for j=1:1:sum(outbound)
                notvalid = 1;
                while notvalid
                    temp_para = sqrt( b*a^2./sum( (randn(b,1)).^2 )' );
                    
                    if ((bounds(i,1) < temp_para)&&(temp_para < bounds(i,2)))
                        parasim(outboundset(j),2) = temp_para;
                        notvalid = 0;
                    end
                end
            end

end