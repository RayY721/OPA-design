function W_result = findbestsolution(W,S,dB,FoV,win)

N = size(W,2);
L = size(S,1);

x = linspace(FoV(1),FoV(2),L)';
ind = find(x >= win(1) & x <= win(2));
S(ind,:) = [];
S_side = S;

% S_side = adjustSmatrix(S,[win.leftend win.rightend],[-1 1],L,0);
max_index = 0;
max_value = zeros(N,1);
for i = 1:1:N
    max_value(i) = max(20*log10(abs(S_side*W(:,i))./max(abs(S*W(:,i)))));
end

max_value(max_value > dB) = -100;
[~,max_index] = max(max_value);

W_result = W(:,max_index);

end