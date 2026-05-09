
function xq = quant(x, B)
% QUANT Uniform B-bit mid-rise quantizer for unity-gain signals
%   xq = quant(x, B) quantizes input signal x (assumed in [-1, 1])
%   using B bits with rounding. Output levels are centered and
%   error is bounded by Q/2.

    % Step size
    Q = 2^(-(B-1));

    % Perform quantization, rounding to nearest level
    xq = Q * (floor(x / Q) + 0.5);

    % Handle edge case at x = 1, keep within valid range
    xq(x >= 1) = 1 - Q/2;
    xq(x <= -1) = -1 + Q/2;
end

