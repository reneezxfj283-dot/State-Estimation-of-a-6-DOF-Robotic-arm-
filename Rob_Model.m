function Mdl = Rob_Model()
rob_Para;
Mdl.trs = trs;

Mdl.Ad{1} = Ad1;
Mdl.Bd{1} = Bd1;
Mdl.Q{1} = Q;
Mdl.I{1} = [200;200;200;55;55;55];

Mdl.Ad{2} = Ad2;
Mdl.Bd{2} = Bd2;
Mdl.Q{2} = Q;
Mdl.I{2} = [200;200;200;55;55;55];

% Mdl.H{1}=[eye(6), zeros(6,6)];
Mdl.H{1}=[1 0 0 0 0 0 0 0 0 0 0 0;     %1
    0 1 0 0 0 0 0 0 0 0 0 0;     %2
    0 0 1 0 0 0 0 0 0 0 0 0;     %3
    0 0 0 1 0 0 0 0 0 0 0 0;     %4
    0 0 0 0 1 0 0 0 0 0 0 0;     %5
    0 0 0 0 0 1 0 0 0 0 0 0 ];    %6
Mdl.R{1} = R;

Mdl.H{2}=[eye(6), zeros(6,6)];
Mdl.R{2} = R;

end