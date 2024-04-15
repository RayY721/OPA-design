function selected_vector = thresholdingw(w, percentage)
    totalenergy = w'*w;
    
    [~, sorted_indices] = sort(abs(w), 'descend');

    accumulated_sum = 0;
    selected_indices = [];
    for i = 1:length(sorted_indices)
        index = sorted_indices(i);
        accumulated_sum = accumulated_sum + w(index)^2;
        selected_indices = [selected_indices, index];
        if accumulated_sum / totalenergy >= percentage
            break;
        end
    end

    selected_vector = zeros(size(w));
    selected_vector(selected_indices) = w(selected_indices);
end