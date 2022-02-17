function label_user = reduced_cluster(B, label_bs)
numU = size(B,2);
label_user = zeros(numU,1);
    for i = 1:numU
        connect_status = B(:,i);
        involved_bs = find(connect_status > 0);
%         involved_bs = np.array([i for i,v in enumerate(connect_status) if v != 0])

        if isempty(involved_bs)
            % the user doesn't connect to any base station
            continue;
        else
            index_cluster = label_bs(involved_bs);
            fre_of_cluster = unique(index_cluster);
            
            if length(fre_of_cluster) == 1
                label_user(i) = index_cluster(1);
                continue;
            else
                clu_weight = zeros(length(fre_of_cluster),1);
                for clu = 1:length(fre_of_cluster)
                    
                    temp = fre_of_cluster(clu);

%                     indicator = [i for i,v in enumerate(index_cluster) if v == temp]
                    indicator = find(index_cluster == temp);
                    base_in_same_clu = involved_bs(indicator);
                    weight = 0;
                    for j =1:length(base_in_same_clu)
                        weight = weight + B(base_in_same_clu(j),i);
                    end
                    clu_weight(clu) = weight;
                end
                max_weight = max(clu_weight);
%                 idx = [i for i,v in enumerate(clu_weight) if v == max_weight]
                idx = find(clu_weight == max_weight);
                if length(idx) == 1
                    label_user(i) = fre_of_cluster(idx);
                else
                    label_user(i) = fre_of_cluster(randi(max(idx)));
                end
            end
        end
    end
end

