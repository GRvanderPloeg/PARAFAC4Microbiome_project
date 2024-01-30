function inputDataPlot_synth(df, timepoints, low_responders, mid_responders, high_responders, speciesName1, speciesName2, speciesName3, overallTitle)
I = size(df, 1);
J = size(df, 2);
K = size(df, 3);

subplot(3,1,1); hold on;
plot(timepoints,reshape(df(low_responders,1,:),size(low_responders,1),K)',"-", "color", [0 1 0]);
plot(timepoints,reshape(df(mid_responders,1,:),size(mid_responders,1),K)',"-", "color", [0 0 1]);
plot(timepoints,reshape(df(high_responders,1,:),size(high_responders,1),K)',"-", "color", [1 0 0]);
title(speciesName1);

% Legend
h = zeros(3, 1);
h(1) = plot(NaN,NaN,'s', "MarkerEdgeColor", [0 1 0], "MarkerFaceColor", [0 1 0]);
h(2) = plot(NaN,NaN,'s', "MarkerEdgeColor", [0 0 1], "MarkerFaceColor", [0 0 1]);
h(3) = plot(NaN,NaN,'s', "MarkerEdgeColor", [1 0 0], "MarkerFaceColor", [1 0 0]);
labels = ["Low responders", "Mid responders", "High responders"];
legend(h, labels, "Location", 'southoutside', "Orientation", "horizontal");

subplot(3,1,2); hold on;
plot(timepoints,reshape(df(low_responders,2,:),size(low_responders,1),K)',"-", "color", [0 1 0]);
plot(timepoints,reshape(df(mid_responders,2,:),size(mid_responders,1),K)',"-", "color", [0 0 1]);
plot(timepoints,reshape(df(high_responders,2,:),size(high_responders,1),K)',"-", "color", [1 0 0]);
title(speciesName2);

% Legend
h = zeros(3, 1);
h(1) = plot(NaN,NaN,'s', "MarkerEdgeColor", [0 1 0], "MarkerFaceColor", [0 1 0]);
h(2) = plot(NaN,NaN,'s', "MarkerEdgeColor", [0 0 1], "MarkerFaceColor", [0 0 1]);
h(3) = plot(NaN,NaN,'s', "MarkerEdgeColor", [1 0 0], "MarkerFaceColor", [1 0 0]);
labels = ["Low responders", "Mid responders", "High responders"];
legend(h, labels, "Location", 'southoutside', "Orientation", "horizontal");

subplot(3,1,3); hold on;
plot(timepoints,reshape(df(low_responders,3,:),size(low_responders,1),K)',"-", "color", [0 1 0]);
plot(timepoints,reshape(df(mid_responders,3,:),size(mid_responders,1),K)',"-", "color", [0 0 1]);
plot(timepoints,reshape(df(high_responders,3,:),size(high_responders,1),K)',"-", "color", [1 0 0]);
title(speciesName3);

% Legend
h = zeros(3, 1);
h(1) = plot(NaN,NaN,'s', "MarkerEdgeColor", [0 1 0], "MarkerFaceColor", [0 1 0]);
h(2) = plot(NaN,NaN,'s', "MarkerEdgeColor", [0 0 1], "MarkerFaceColor", [0 0 1]);
h(3) = plot(NaN,NaN,'s', "MarkerEdgeColor", [1 0 0], "MarkerFaceColor", [1 0 0]);
labels = ["Low responders", "Mid responders", "High responders"];
legend(h, labels, "Location", 'southoutside', "Orientation", "horizontal");

sgtitle(overallTitle);