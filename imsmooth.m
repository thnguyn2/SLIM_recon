%smoothing images by n pixels

function A = imsmooth(M, n)

A = conv2(M, ones(n), 'same')/n^2;
