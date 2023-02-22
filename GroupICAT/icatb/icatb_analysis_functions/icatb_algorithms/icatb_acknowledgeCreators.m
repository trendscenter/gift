function icatb_acknowledgeCreators
% icatb_acknowledgeCreators

icatb_defaults;

global FLAG_ACKNOWLEDGE_CREATORS;

if strcmp(lower(FLAG_ACKNOWLEDGE_CREATORS), 'on')
    %setup strings
    introductionString = ['The ICA algorithms used in this toolbox were created by multiple individuals. Please acknowledge the creators of the algorithms.'];
    infomaxString = str2mat('', '1. Infomax Algorithm: Information Maximization', 'Created By: Scott Makeig', '');

    fastString1 = '2. Fast ICA Algorithm:';
    fastString2 = 'Created By: Hugo Gavert, Jarmo Hurri, Jaakko Sarel & Aapo Hyvarinen.';
    fastString = str2mat(fastString1, fastString2, '');

    ericaString1 = '3. ERICA Algorithm: Equivariant Robust Indepedent Component Analysis algorithm';
    ericaString2 = 'Created By: Sergio Cruces, Luis Castedo & Andrzej Cichocki.';
    ericaString = str2mat(ericaString1, ericaString2, '');


    simbecString1 = '4. Simbec Algorithm: Simultaneous Blind signal Extraction using Cumulants';
    simbecString2 = 'Created By: Sergio Cruces, Andrzej Cichocki & S. Amari.';
    simbecString = str2mat(simbecString1, simbecString2, '');


    evdString1 = '5. EVD: EVD method for source separation';
    evdString2 = 'Created By: Andrzej Cichocki & Pando Georgiev.';
    evdString = str2mat(evdString1, evdString2, '');

    jadeopString1 = '6. JADEOP: Jadeop Algorithm:';
    jadeopString2 = 'Created By: Cardoso and Optimized by: JY. Terazono & A. Cichocki.';
    jadeopString = str2mat(jadeopString1, jadeopString2, '');

    amuseString1 = '7. Amuse Algorithm: BSS using singular value decomposition';
    amuseString2 = 'Created By: A. Cichocki.';
    amuseString = str2mat(amuseString1, amuseString2, '');

    optimalString1 = '8. SDD ICA: Source Density Driving ICA';
    optimalString2 = 'Created By: Baoming Hong & Vince Calhoun.';
    optimalString = str2mat(optimalString1, optimalString2, '');

    semiblindString1 = '9. SBICA: Semi-blind Infomax ICA Algorithm';
    semiblindString2 = 'Created By: Vince Calhoun.';
    semiblindString = str2mat(semiblindString1, semiblindString2);

    creditString = str2mat(introductionString, infomaxString, fastString, ericaString, simbecString, ...
        evdString, jadeopString, amuseString, optimalString, semiblindString);

    for i = 1:size(creditString, 1)
        msgString{i} = creditString(i, :);
    end

    %creditString = str2mat(introductionString, infomaxString, ericaString, optimalString, fastString);

    % make user acknowledge algorithm creators
    % answer = questdlg(creditString,'ICA Algorithm Credit','Acknowledge','More Info','More Info');
    %
    % if(strcmp(answer,'More Info'))
    %     disp('The More Info feature will be added later');
    % end

    %CreateMode.WindowStyle = 'modal'; CreateMode.Interpreter = 'none';
    %msg_handle = msgbox(creditString, 'ICA Algorithm Credit', 'help', CreateMode);

    msg_handle = icatb_dialogBox('title', 'ICA Algorithm Credit', 'textBody', msgString, 'textType', 'large');

    waitfor(msg_handle);

end