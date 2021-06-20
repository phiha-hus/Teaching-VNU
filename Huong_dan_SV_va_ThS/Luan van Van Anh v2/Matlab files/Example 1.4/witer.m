function	eps	=	witer(W)
global	h	A	Ad
eps = W*expm(W+h*A)-h*Ad;
end