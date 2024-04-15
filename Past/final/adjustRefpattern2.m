function F_new = adjustRefpattern2(F,FoV,win,L,loose_num)
x = linspace(FoV(1),FoV(2),L)';
ind = find((0 < win(1) - x) & (win(1) - x <= loose_num) | (0 < x - win(2)) & (x - win(2) <= loose_num));
F(ind,:) = [];
F_new = F;
end