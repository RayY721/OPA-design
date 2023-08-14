function S_new = adjustSmatrix(S,left,right,win,L,loose_num)
% loose_num in degree - The range that will be neglected

x = linspace(left,right,L)';
ind = find((abs(x - win(1)) <= loose_num) | (abs(x - win(2)) <= loose_num));
S(ind,:) = [];
S_new = S;
end