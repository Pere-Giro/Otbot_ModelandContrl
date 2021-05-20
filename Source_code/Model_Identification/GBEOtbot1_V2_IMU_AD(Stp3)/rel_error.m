function [rel_error,abs_error] = rel_error(real_val,measured_val)
%REL_ERROR Summary of this function goes here
%   This function computes and returns the value of relative error computed
%   as abs(real value - measured value)/ real value * 100. If the real
%   value is 0 then it returns the absolute error without multiplying
lsv = max(size(real_val));

rel_error = zeros(lsv,1);
abs_error = zeros(lsv,1);

for i=1:lsv
    if real_val(i,1)~=0
        rel_error(i,1) = abs(real_val(i,1) - measured_val(i,1))./abs(real_val(i,1))*100;
    else
        rel_error(i,1) = abs(real_val(i,1) - measured_val(i,1));
    end
    abs_error(i,1) = abs(real_val(i,1) - measured_val(i,1));
end
end

