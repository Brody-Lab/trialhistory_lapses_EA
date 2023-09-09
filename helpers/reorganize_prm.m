function[d] = reorganize_prm(rats, prm_field)

type = fields(rats(1).(prm_field));
prm_type = fields(rats(1).(prm_field).(type{1}));

for r = 1:length(rats)
    for t = 1:length(type)
        for p = 1:length(prm_type)
            d.(type{t}).(prm_type{p})(r) = rats(r).(prm_field).(type{t}).(prm_type{p});
        end
    end
    d.ntrials(r) = rats(r).ntrials;
    d.accuracy(r) = rats(r).accuracy;
end
end