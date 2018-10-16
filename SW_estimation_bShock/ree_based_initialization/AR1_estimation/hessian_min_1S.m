% Copyright (C) 2001 Michel Juillard
% Copyright (C) 2007 Sergey Slobodyan
%
% computes 1-sided second order partial derivatives
% uses Abramowitz and Stegun (1965) formulas 25.3.24 and 25.3.27 p. 884
function [hessian_mat,info,x_min,f_min] = hessian_min(func,x,scale)
lik_penalty=0;
gstep_=0.01;
cliff_step = max(0.5,0.9 * lik_penalty / 2);

scale_ = scale;
%func = str2func(func);
n=size(x,1);
%h1 = scale_ * max(abs(x),sqrt(gstep_)*ones(n,1))*eps^(1/6);
h1=max(abs(x),gstep_*ones(n,1))*eps^(1/3);
done = 0;
x_min = x;
f0 = 1e+8;
f_min = 1e+8;
while (done == 0)
    cliff = 0;
    while (cliff == 0)
        cliff_loc = zeros(1,n);
        if f0 > f_min
            x = x_min;
        end
        h_1=h1;
        xh1=x+h1;
        h1=xh1-x;
        xh1=x-h_1;
        h_1=x-xh1;
        xh1=x;
        info = 0;

%         f0=feval(func,x,varargin{:});
f0=feval(func,x);
        f_min = f0;
        f1=zeros(size(f0,1),n);
        f_1=f1;
        hessian_mat = zeros(n,n);

        for i=1:n
            xh1(i)=x(i)+h1(i);
%             f1(:,i)=feval(func,xh1,varargin{:});
 f1(:,i)=feval(func,xh1);
            xh1(i)=x(i)-h_1(i);
%             f_1(:,i)=feval(func,xh1,varargin{:});
 f_1(:,i)=feval(func,xh1);
            xh1(i)=x(i);
            i=i+1;
        end

        % forming the vector of 'good', or cliff-free, directions
        [trash,index] = min([f1;f_1]);
        hk = zeros(1,length(h1));
        hk(find(index==1)) = h1(find(index==1));
        hk(find(index==2)) = -h_1(find(index==2));
        fk = zeros(1,length(h1));
        fk(find(index==1)) = f1(find(index==1));
        fk(find(index==2)) = f_1(find(index==2));

        [f_min,ind] = min(fk);
        if f_min < f0
            x_min = x;
            x_min(ind) = x(ind) + hk(ind);
            info = 1;
        end

        temp = (fk - f0) > cliff_step;
        if ~isempty(find(temp, 1))
            h1(find(temp)) = h1(find(temp)) / 10;
            cliff = 1;
            break
        end

        for i=1:n
            if i > 1
                hessian_mat(i,1:i-1) = hessian_mat(1:i-1,i);
            end
            xh1 = x;
            xh1(i) = x(i) + 2 * hk(i);
%             t1 = feval(func,xh1,varargin{:});
 t1 = feval(func,xh1);
            if (t1 - f0) > cliff_step
                cliff_loc(i) = 1;
                cliff = 1;
            end
            if t1 < f_min
                f_min = t1;
                x_min = xh1;
                info = 1;
            end
            hessian_mat(i,i)=(f0 + t1 -2*fk(i))./(hk(i)^2);
            for j=i+1:n
                xh1(i) = x(i) + hk(i);
                xh1(j) = x(j) + hk(j);
%                 t1 = feval(func,xh1,varargin{:});
 t1 = feval(func,xh1);
                if (t1 - f0) > cliff_step
                    cliff_loc(i) = 1;
                    cliff_loc(j) = 1;
                    cliff = 1;
                end
                if t1 < f_min
                    f_min = t1;
                    x_min = xh1;
                    info = 1;
                end
                hessian_mat(i,j) = ( t1 + f0 - fk(i) - fk(j) )./(hk(i)*hk(j));
                xh1(i)=x(i);
                xh1(j)=x(j);
            end
        end
        cliff_loc
        h1(find(cliff_loc)) = h1(find(cliff_loc)) / 10;
        eig(hessian_mat)' %#ok<NOPRT>
        if cliff == 1 % reduce the step
            scale_ = scale_ / 2;
            if scale_ / scale < 1.e-2 % the step is too low
                done = 1;
                info = 2;
            end
        elseif min(eig(hessian_mat)) < 0 % Hessian is not positively definite
            [V,D] = eig(hessian_mat);
            [trash,ind] = min(diag(D));
            [alfa_min,f_min_] = fmincon(@Fun_min,0,[],[],[],[],-1,1,[],[]);
            if f_min_ < f_min
                f_min = f_min_;
                x_min = x + alfa_min * max(xh1) * V(:,ind);
            end
            done = 1;
            info = 1;
            scale_ = scale_ / 2;
            if scale_ / scale < 1.e-2 || f0 - f_min < 1.e-06 % the step is too low or improvement too small
                done = 1;
                info = 2;
            end            
            alfa_min
            f_min
            trash
        else
            done = 1;
            info = 1;
        end
        break
    end
end

    function f = Fun_min(alfa)
%         f = feval(func,x + alfa*max(xh1)*V(:,ind),varargin{:});
  f = feval(func,x + alfa*max(xh1)*V(:,ind));
    end

end


% 11/25/03 SA Created from Hessian_sparse (removed sparse)
% 05/11/2007 Sergey Modified to accomodate bad likelihood structure of models with
% adaptive learning
