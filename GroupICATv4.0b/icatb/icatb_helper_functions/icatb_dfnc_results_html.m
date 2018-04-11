function icatb_dfnc_results_html(dfncInfo)
%% Print dfnc results to HTML file
%

helpStr = 'Creating HTML file. This will involve writing jpeg files to the disk. Please wait ...';
helpH = helpdlg(helpStr);
disp(helpStr);

html_file = fullfile(dfncInfo.outputDir, 'html', [dfncInfo.prefix, '_results.html']);
start_string = '<html><head><title> DFNC Results </title></head>';
start_string = [start_string, '<p> </p><h2> Contents </h2><p><ul>'];

titles = {'FNC Oscillations','Cluster Stats', 'Meta State Analysis'};

for nR = 1:length(titles)
    start_string = [start_string, '<li><h3><a href="#results', num2str(nR), '">', titles{nR}, '</a></h3></li>'];
end

start_string = [start_string, '</ul></p><p> </p>'];

results = icatb_dfnc_results(dfncInfo, 'fnc oscillations');

titleStr = titles{1};

results_string1 = ['<hr> <h2 align = "center">', '<a name="results1">', titleStr, '</a></h2>'];
for nR = 1:length(results);
    results_string1 = [results_string1, '<p align = "center">', results(nR).text, '</p>'];
    results_string1 = [results_string1, '<p align = "center"> <img src = "', results(nR).file, '" > </img>'];
end

end_string =  '</html>';

results = icatb_dfnc_results(dfncInfo, 'clusters');

titleStr = titles{2};

results_string2 = ['<hr> <h2 align = "center">', '<a name="results2">', titleStr, '</a></h2>'];
for nR = 1:length(results);
    results_string2 = [results_string2, '<p align = "center">', results(nR).text, '</p>'];
    results_string2 = [results_string2, '<p align = "center"> <img src = "', results(nR).file, '" > </img>'];
end


try
    results = icatb_dfnc_results(dfncInfo, 'meta state analysis');
    titleStr = titles{3};
    
    results_string3 = ['<hr><h2 align = "center">', '<a name="results3">', titleStr, '</a></h2>'];
    for nR = 1:length(results);
        results_string3 = [results_string3, '<p align = "center">', results(nR).text, '</p>'];
        results_string3 = [results_string3, '<p align = "center"> <img src = "', results(nR).file, '" > </img>'];
    end
    results_string = [start_string,  results_string1, results_string2, results_string3, end_string];
catch
    results_string = [start_string,  results_string1, results_string2, end_string];
end


%results_string = [start_string,  results_string1, results_string2, end_string];

dlmwrite(html_file, results_string, '');

icatb_openHTMLHelpFile(html_file);

try
    delete(helpH);
catch
end

fprintf('Done\n');
