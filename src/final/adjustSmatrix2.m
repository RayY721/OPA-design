function S_new = adjustSmatrix2(S,FoV,win,L,loose_num)
% loose_num in degree - The range that will be neglected

x = linspace(FoV(1),FoV(2),L)';
ind = find((0 < win(1) - x) & (win(1) - x <= loose_num) | (0 < x - win(2)) & (x - win(2) <= loose_num));
S(ind,:) = [];
S_new = S;
end