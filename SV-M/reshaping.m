%define a function to reshape the matrices after the for loop
%initially, each row contains all the quantiles for 4 successive weeks
%we now split each row into 4 separate rows that corresponds to the 4 weeks
function reshaped_mat = reshaping(mat,iter)

reshaped_mat = zeros(4*(iter+1),6);
for i = 0:iter
    for j = 0:3
        start_point = 6*j+1;
        end_point = 6*(j+1);
        reshaped_mat(4*i+j+1,:) = mat(i+1,start_point:end_point);
    end
end

end