%function [x,fval]=mylinprog(c,A,b,Aeq,beq,vlb,vub)
%[x,fval]=linprog(c,A,b,Aeq,beq,vlb,vub);
%end
function [x,fval]=mylinprog(c,A,b,vlb,vub)
[x,fval]=linprog(c,A,b,[],[],vlb,vub);
end
