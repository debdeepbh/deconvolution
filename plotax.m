[ax, aximp] = prepsig(5);
%plot(ax')
tot = 5
for i = 1:tot
    subplot(tot,2, 2*i-1)
    plot(ax'(:,i))
    axis('tight')
    axis('nolabel')

    subplot(tot,2, 2*i)
    plot(aximp'(:,i))
    %axis('off')
    axis('tight')
    axis('nolabel')
end

%demo_save(1, 'anitadata', '800,600')

