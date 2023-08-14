function F_new = adjustrefpattern(F,left,right,win,L,loose_num)
x = linspace(left,right,L)';
ind = find((abs(x - win(1)) <= loose_num) | (abs(x - win(2)) <= loose_num));
F(ind,:) = [];
F_new = F;
end