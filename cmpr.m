function ix=cmpr(c1,c2)
% compare c1 with every elemnt of c2
% c1 and c2 are cell
L=length(c2);
[i0,i1]=size(c1);
ix=zeros(i0,i1);
for l=1:L
    ix= ix + strcmp(c1,c2(l));
end
end