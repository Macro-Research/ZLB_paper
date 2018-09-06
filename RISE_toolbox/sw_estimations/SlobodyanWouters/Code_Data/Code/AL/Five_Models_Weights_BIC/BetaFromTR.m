function [beta,SecMom,R_beta] = BetaFromTR(T,R,Q)

%  This program was written for Sergey Slobodyan and Raf Wouters DSGE
%  estimation under adaptive learning project.
%  It takes transmission mechanism of the model (typically, derived at some 
%  REE) and forms agents' beliefs which are consistent with the stationary
%  distribution implied by this transmission mechanism
%  Sergey.Slobodyan@cerge-ei.cz, August 18, 2011

global options_ 

crit = 1.e-08;

ys_list = options_.ys_list;
shock_list = options_.shock_list;
yf_list = options_.yf_list;

%%%% Deriving the theoretical  betamat for constant gain learning. 
%%%% Rows are arranged alphabetically. 
%%%% Columns are ENDOGENOUS vars first, SCHOCKS second, alpha within groups. 

tmp = lyapunov_symm(T,R*Q*R');
tmp = 0.5 * (tmp + tmp');

XX11 = tmp(ys_list,ys_list);
XX22 = tmp(shock_list,shock_list);
XX21 = T(shock_list,:)*tmp(:,ys_list);
XX = [XX11 XX21'; XX21 XX22];
XX = 0.5 * (XX + XX');

% Re-scaling
C = diag(1 ./ sqrt(diag(XX)) );
SecMom = XX;
XX = C * XX * C;

XY1 = T(yf_list,:) * tmp(:,ys_list);
XY2 = tmp(yf_list,shock_list);
XY = [XY1 XY2];

beta = C * (XX \ (XY * C)');

options_.tmp = tmp;

%%%% Beliefs and \Sigma matrix for the Kalman filter learning

if options_.kalman_algo > 500
    for nm = 1:options_.num_mod
        %%%% So far, it's only for the states!
        ys_lists = options_.m(nm).y_st_full;
        ys_list_all = nonzeros(ys_lists);
        Mask = zeros(length(ys_list_all),size(ys_lists,2));
        Mask_const = zeros(length(ys_list_all)+size(ys_lists,2),size(ys_lists,2));
        offset = 0;
        offset_const = 0;
        for i = 1:size(ys_lists,2)
            len = length(nonzeros(ys_lists(:,i)));
            Mask(offset+1:offset+len,i) = 1;
            Mask_const(offset_const+1:offset_const+len+1,i) = 1; %#ok<AGROW>
            offset = offset + len;
            offset_const = offset_const + len + 1;            
        end
        if options_.kalman_algo > 600
            options_.m(nm).Mask_const = Mask_const;
        end
        XX_MA = tmp(ys_list_all,ys_list_all) .* (Mask * Mask');
        XX_MA = 0.5 * (XX_MA + XX_MA');
        XY_MA = T(yf_list,:) * tmp(:,ys_list_all);
        XY_MA_vec= Mask .* XY_MA';
        XY_MA_vec = nonzeros(XY_MA_vec(:));

        beta = XX_MA \ XY_MA_vec;
        XY_MA_ = kron(Mask',ones(length(yf_list),1));
        offset = 0;
        for i = 1:length(yf_list)
            XY_MA_(offset+1:offset+length(yf_list),:) = XY_MA .* XY_MA_(offset+1:offset+length(yf_list),:);
            offset = offset+length(yf_list);
        end

        beta_MA = kron(beta,ones(1,length(yf_list))) .* Mask;
        R_beta = tmp(yf_list,yf_list) - reshape(XY_MA_ * beta,length(yf_list),length(yf_list))...
            -transpose(reshape(XY_MA_ * beta,length(yf_list),length(yf_list)))...
            + beta_MA' * tmp(ys_list_all,ys_list_all) * beta_MA;
        R_beta = 0.5 * (R_beta + R_beta');
        %%%% Creating a positively definite matrix
        [V,D] = eig(R_beta);
        DD = diag(D);
        if length(find(DD < crit))
            DD(DD < crit) = min(DD(DD > crit));
            R_beta = V * diag(DD) / V;
            R_beta = 0.5 * (R_beta + R_beta');
        end
        inv_V = Mask * (R_beta \ Mask');
        options_.m(nm).Q_beta = inv(tmp(ys_list_all,ys_list_all) .* inv_V);
        options_.m(nm).Mask = Mask;
        options_.m(nm).beta = beta;
        options_.m(nm).R_beta = R_beta;

    end
elseif options_.kalman_algo > 300
%%%% Deriving R_beta. Currently I assume that it's MSV and so I have to
%%%% change the timing assumptions

    %%%% Making shocks appear as t-1 vars, not t vars
    ys_list = [ys_list; shock_list];
    shock_list = [];

    XX11 = tmp(ys_list,ys_list);
    XX22 = tmp(shock_list,shock_list);
    XX21 = T(shock_list,:)*tmp(:,ys_list);
    XX_ = [XX11 XX21'; XX21 XX22];
    XX_ = 0.5 * (XX_ + XX_');

    % Re-scaling
    C = diag(1 ./ sqrt(diag(XX_)) );
    XX_ = C * XX_ * C;

    XY1 = T(yf_list,:) * tmp(:,ys_list);
    XY2 = tmp(yf_list,shock_list);
    XY_ = [XY1 XY2];
    XY_ = XY_ * C;

    %%%% Given that the matrix T is singular for the MSV, have to use pinv
    beta_ = C * pinv(XX_) * XY_';

    R_beta = tmp(yf_list,yf_list) - beta_' * XY';
    R_beta = 0.5 * (R_beta + R_beta');
    %%%% Creating a positively definite matrix
    [V,D] = eig(R_beta);
    DD = diag(D);
    DD(DD < crit) = min(DD(DD > crit));
    R_beta = V * diag(DD) / V;
    R_beta = 0.5 * (R_beta + R_beta');
else
    R_beta = [];
end
