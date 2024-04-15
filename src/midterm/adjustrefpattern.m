function F_new = adjustrefpattern(F,FoV,win,L,loose_num)
x = linspace(FoV(1),FoV(2),L)';
ind = find((abs(x - win(1)) <= loose_num) | (abs(x - win(2)) <= loose_num));
F(ind,:) = [];
F_new = F;
end