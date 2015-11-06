function C = sense_lin(t,x,S)
g=S.g;

syms xsym;
gsym=@(xsym)g(xsym);

gd=diff(gsym(xsym));
gfun = matlabFunction(gd);
try
  gdash_val=gfun(x);
catch ME
  if (strcmp(ME.identifier ,'MATLAB:TooManyInputs'))     
   gdash_val=gfun();
  end
end

C=[-gdash_val*S.omega*sin(S.omega*t) g(cos(S.omega*t))];
end