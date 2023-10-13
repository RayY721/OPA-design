function S_new = adjustSmatrix(S,FoV,win,L,loose_num)
% loose_num in degree - The range that will be neglected

x = linspace(FoV(1),FoV(2),L)';
ind = find((abs(x - win(1)) <= loose_num) | (abs(x - win(2)) <= loose_num));
S(ind,:) = [];
S_new = S;
end