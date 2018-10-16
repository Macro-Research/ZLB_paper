            try
%largest_eig2(tt)=pp_filtered(tt)*(eigs(gamma1_1,1))+(1-pp_filtered(tt))*(eigs(gamma1_2,1));
largest_eig2(tt)=max(abs(eigs(gamma1_1,1)),abs(eigs(gamma1_2,1)));
%largest_eig2(tt)=((eigs(gamma1_2,1)));
            catch
                
                largest_eig2(tt)=1.01;
            end