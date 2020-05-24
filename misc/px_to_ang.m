function ret = px_to_ang(simdim,kev)
if ~exist('kev','var')
      kev = 300;
end
lambda = 12.3986./sqrt((2*511.0+kev).*kev) * 10^-10;
ret = lambda*1/(2*simdim/1000)*10^10; 

end