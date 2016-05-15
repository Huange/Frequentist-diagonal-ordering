function triu_new_adj_matrix = frequentist_diagonal_ordering(orig_matrix,n_iter)
% Input attributes: orig_matrix - the original adjacent matrix before frequentist diagonal ordering
%                   n_iter - the predefined number of iterations
% Output attributes: triu_new_adj_matrix - upper triangular of the new adjacent matrix
% key intermediate attributes: OUTPUT_order - a "1"x"n" cell array; each element represents an updated marker order
%                              OUTPUT_matrix - a "1"x"n" cell array; each element represents an updated adjacent matrix
%                              OUTPUT_score - a "1"x"n" cell array; each element represents the total number of "1" entries 
%                                             on the first super diagonal of an updated adjacent matrix
% Author: Huange Wang
% Create date: 2016-05-15

n = size(orig_matrix,2);
OUTPUT_score = cell(1,n_iter);
OUTPUT_order = cell(1,n_iter);
OUTPUT_matrix = cell(1,n_iter);

for N = 1:n_iter
    
    matrix = orig_matrix;
    L = fliplr(1:n);
    
    orig_order = 1:n;
    output_order = zeros(0);
    left_order = orig_order(~ismember(orig_order,output_order));
    
    output_matrix = zeros(n);
    output_score = 0;
    
    for i = 1:n
        l = L(i);
        test_score = zeros(1,l);
        test_column = zeros(n,l);
        
        for t = 1:l
            A = matrix(:,t);
            
            if i > 1
               for j = 1:(i-1)
                   test_column(j,t) = A(output_order(j));
               end
            end
            B = A(~ismember(orig_order,output_order));
            test_column(i,t) = B(t);
            
            if i > 1
               test_score(t) = test_column(i-1,t);
            end
        end

        max_score = max(test_score);
        output_score = output_score + max_score;
        
        short_list = left_order(ismember(test_score, max_score));
        if length(short_list) > 1
           select = randsample(short_list,1);
        else
           select = short_list;
        end
        output_order(i) = select;
        
        K = 1:l;
        k = K(ismember(left_order,output_order(i)));
        matrix = matrix(:,K(~ismember(K,k)));    
        output_matrix(:,i) = test_column(:,k);
        left_order = orig_order(~ismember(orig_order,output_order));  
        
        OUTPUT_score{N} = output_score;
        OUTPUT_order{N} = output_order;
        OUTPUT_matrix{N} = output_matrix;
    end
end

M = max([OUTPUT_score{1,:}]);
list = find([OUTPUT_score{1,:}] == M);
n_list = length(list);
list_order = cell(n_list,1);
for i = 1:n_list
    list_order{i} = OUTPUT_order{list(i)};
end

n = size(OUTPUT_order{1,1},2);
neighbors = cell(n,1);

for t = 1:n_list
    list = list_order{t};
    for i = 1:n
        mark = find(list == i);
        L = length(neighbors{i});
        if  mark == 1
            neighbors{i}(1,L+1) = list(mark+1);
        elseif mark == n
            neighbors{i}(1,L+1) = list(mark-1);
        elseif mark > 1 && mark < n
            neighbors{i}(1,L+1) = list(mark-1);
            neighbors{i}(1,L+2) = list(mark+1);
        end
    end
end

freq = zeros(n);
for i = 1:n
    for j = 1:n
        freq(i,j) = length(find(neighbors{i} == j));
    end
end
freq = freq/n_list;          % freq is symmetric !

neighbors_short = cell(n,1);
for i = 1:n
    for j = 1:n
        L = length(neighbors_short{i});
        if freq(i,j) >= 0.5               % the threshold of 0.5 is adjustable!
           L = L + 1;
           neighbors_short{i}(1,L) = j;
        end
    end
end

new_adj_matrix = zeros(n);
for i = 1:n;
    x = neighbors_short{i,1};
    L = length(x);
    if L > 0
       for l = 1:L
           new_adj_matrix(i,x(l)) = 1;
       end
    end
end
triu_new_adj_matrix = triu(new_adj_matrix);
