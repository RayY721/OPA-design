%%
function reference_pattern = desiredbeam(left,right,win,L)   % only symmetry rect 
% sll in dB, 
% left and right indicate the range of the main lobe
% L is the sampling number with the whole FOV
    % sll_value = 10^(sll/20);
    
    syms y(x)
    y(x) = piecewise((x >= win(1) & x <= win(2)),1,0);  
    x = linspace(left,right,L)';
    reference_pattern = double(y(x));
end