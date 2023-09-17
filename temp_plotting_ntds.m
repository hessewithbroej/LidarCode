%% plot generation. Run main first and do not clear all


n_tds = [0, 10.^(-2:1:2)];

outputs = arrayfun(@(ntd) temp_function(lowest_constituent,ntd,par_detected_binned), n_tds, 'UniformOutput',false);

%%

colors = [1,0,0;
    1,0.5,0;
    1,1,0;
    0.5,1,0;
    0,1,0;
    0,1,0.5;
    0,1,1;
    0,0.5,1;
    0,0,1;
    0.5,0,1;
    1,0,1;
    1,0,0.5;];
colororder(colors)



figure
hold on
for i =1:numel(outputs)
    plot(0:.001:1, fliplr(outputs{i}), 'color',colors(i,:))
end

xlabel("$P_{n,receive}$", 'Interpreter', 'latex')
ylabel("$P_{n,detect}$", "Interpreter", "latex")
title("$P_{n,receive}$ vs $P_{n,detect}$ for various $n_{td}$", "Interpreter", "latex")
legend(string(n_tds), 'Location','northwest')

