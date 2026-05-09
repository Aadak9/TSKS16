function y = AGCunity(x)
% Normalize signal to range [-1, 1]

y = x / max(abs(x));

end