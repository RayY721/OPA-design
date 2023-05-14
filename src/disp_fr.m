function disp_fr(w,N,numG)
for i = 1:1:numG
    X = ['The filling rate of No.', num2str(i), ' is ', num2str(sum(w(:,i))/(N/2)*100)];
    disp(X)
end
end