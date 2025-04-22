function[M_est,M_error] = curve_fit(x,y,err)

%make sure the input data is a column vector
temp = size(x);
if (temp(1) < temp(2))
    x = x';
end
temp = size(y);
if (temp(1) ==1)
    y = y';
end


%if no error vector is given then use an identity matrix
if (nargin == 2)
    W = diag(ones(1,length(x)));
else
    if length(err) == 1
        err = err*ones(size(y));
    end
    w = 1./(err.^2);
    W = diag(w);
end    


G = [x];

temp = inv(G'*W*G)* G' * W;
M_est =  temp * y;

if (nargin == 2)
    resid = y - G*M_est;
    err = cov(resid) * ((length(x)-1)/(max(size(G))-min(size(G)))) * inv(G' * G);
else
    w = err.^2;
    cov_y = diag(w);
    err = temp * cov_y * temp';
end
M_error = sqrt(diag(err));