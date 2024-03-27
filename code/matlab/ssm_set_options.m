function option = ssm_set_options()

    option.group = false;
    option.filepath = '/Users/jolandamalamud/Jolanda/PhD/ssm';
    option.foldername = 'leemput_patient';
    if option.group 
        option.filename = 'twindata';
    else
        option.filename = 'leemput_patient'; %'leemput_patient' or 'leemput_control'
    end
    option.Adiag = false;
    option.input = false;
    option.bias = false;
    option.pca = false;
    option.gap = 170;
    option.empirical = true;
  
end