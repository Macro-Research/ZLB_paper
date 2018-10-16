            try
%largest_eig2(tt)=max(abs(eigs(gamma1_1,1)),abs(eigs(gamma1_2,1)));
largest_eig2(tt)=(abs(eigs(gamma1_1,1)));

            catch
                largest_eig2(tt)=1.01;
            end