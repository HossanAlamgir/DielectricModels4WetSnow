function styleSet=plot_styles()

% Colorblind-friendly RGB triplets
cbColors = {
    [0, 0, 0],              % Black (neutral and high contrast)
    [0, 114, 178]/255,      % Blue (distinct and common)
    [230, 159, 0]/255,      % Orange (bright and warm)
    [0, 158, 115]/255,      % Bluish Green (cool and distinct)
    [240, 228, 66]/255,     % Yellow (bright, but low contrast on white)
    [86, 180, 233]/255,     % Sky Blue (lighter than main blue)
    [213, 94, 0]/255,       % Vermillion (red-orange, warm)
    [204, 121, 167]/255     % Reddish Purple (more subtle)
};


% Basic line styles
lineStyles = {'-', '--', ':', '-.'};

% Create cell array of structs for plotting
styleSet = {};
for i = 1:length(cbColors)
    for j = 1:length(lineStyles)
        styleSet{end+1} = struct('Color', cbColors{i}, 'LineStyle', lineStyles{j}); %#ok<SAGROW>
    end
end
end