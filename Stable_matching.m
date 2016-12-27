function match_all = stable_matching(v,quota)
% return 0 ~ I-1
J_t = size(v,1);
I = size(v,2);
match_all = -1*ones(J_t,1);
% (1) expand I and v according to quota
I_expand = -1*ones(sum(quota),1);
v_expand = zeros(J_t,sum(quota));
l = 1;
for i = 1:I
    for q = 1:quota(i)
        I_expand(l) = i;
        v_expand(:,l) = v(:,i);
        l = l+1;
    end;
end;
% (2) one-to-one matching - match_expand(J_t)
match_expand = -1*ones(J_t,1);
[B,index] = sort(v_expand,2,'descend');
count = sum(match_expand>-1);
not_proposed = index;
while count < min(J_t, sum(quota))
    free = find(match_expand==-1);
    j = free(1);
    free = find(not_proposed(j,:)>-1);
    i = not_proposed(j,free(1));
    not_proposed(j,find(not_proposed(j,:)==i)) = -1;
    v_temp = v_expand(j,i);
    % i is free
    if size(find(match_expand==i),1) == 0
        match_expand(j) = i;
    else
        j_a = find(match_expand==i);
        v_current = v_expand(j_a,i);
        if v_temp > v_current;
            match_expand(j) = i;
            match_expand(j_a) = -1;
        end;
    end;
    count = sum(match_expand>-1);
end;
% (3) collapse to one-to-many matching - match_all(J_t)
for j = 1:J_t
    if match_expand(j) == -1
        match_all(j) = -1;
    else
        match_all(j) = I_expand(match_expand(j))-1;
    end;
end;
