plot(RespU.(Evar{3}).data);
hold on
for i=1:NPer
plot(IRFsU.Simul(:,i))
end
hold off