# ANOVA Implemented Into dFNC 11/07/2025

A 1-way (N-level) ANOVA has now been implemented in the GIFT GUI. You may find the following group-analysis example helpful: (https://github.com/trendscenter/gift/tree/master/GroupICAT/demo/gift_group_analysis).<BR>
<BR>
In order to run this you need to install the latest version (v4.0.6.15) of GIFT found at:<BR>
https://github.com/trendscenter/gift<BR>
Then you need to click the green [<> Code] button and then select the "Download ZIP" from the list<BR>
<BR>
Below are step-by-step instructions:<BR>
<BR>
GETTING STARTED EXAMPLE<BR>
1. Run your ICA:<BR>
First you need to run some type of GIFT ICA. For this demo we used a batch script such as the batch example found at https://github.com/trendscenter/gift/blob/master/GroupICAT/icatb/icatb_batch_files/batch_constrained_ica.m<BR>
2. Launch GIFT (GUI mode):<BR>
Type gift in the MATLAB command window.<BR>
<BR>
CONTINUING WITH ANOVA MODEL.<BR>
3. In GIFT, click the [Dynamic Functional Connectivity] button.<BR>
4. Click the [Temporal dFNC (ICA)] button.<BR>
5. Click the [Setup / Run Analysis] button.<BR>
6. Browse to your parameter file from your ICA run you want to explore with ANOVA dFNC. Alternativelu, you may download demodata from our group-analysis page (above), picking neuromark_ica_parameter_info.mat from the demo_input3neuromark/gift_out directory.<BR>
Demo data link is: https://trends-public-website-fileshare.s3.amazonaws.com/public_website_files/software/gift/data/demo_input3neuromark.zip<BR>
7. Select output folder your ANOVA results will be saved.<BR>
8. Add components (and domains) by clicking the [+] button. Then type in domain name at Enter Network name and then select the components that domain includes (use shift or control to select multiple components for domain). Finally click [Done] button to add the domain and components. You may also add more than one domain if you choose.<BR>
9. After you added your domains that you want you may click the [Run] button and wait about a minute until the processing is done (seen in matlab command window)<BR>
10. Click [Post-processing] button and select your {prefix}_dnfc.mat file and a Post-processing dfnc parameters window pops up. You may leave all the options enabled, but perhaps you want to change the number of k-means clusters (instead of having 6 clusters). Also, at this same screen we recommend to change the [Cluster options], by pressing [Cluster options] in the upper left, getting a "Figure 1: Select K-means options", changing the "Enter maximum number of iterations" from 150 to 10000 (to find the best global maximum). For the same reason you may also increase the "Number of times to repeat the clustering" from 5 to 75. Clock [OK] to get the new settings.<BR>
11. In the Post-processing dfnc parameters window, click the [Done] button and wait between 1 minute to one hour depending on what settings you picked.<BR>
12. Skipping the display, you may jump to the [Stats'] button and select your {prefix}_dfnc.mat file and select your output folder.<BR>
13. Change the design criteria to "1-way, x-level anova" in the dropdown and then a small "ANOVA Levels" window pops up where you enter number of levels in the textbox, perhaps 3 and then hit the [OK] button (not all ANOVA features work for more than 5 levels).<BR>
14. Enter name for Group 1 (level 1) e.g., Male1, and then select the subjects (using the mouse and the control button on the keyboard) to select subjects to add. You may also enter the comma separated vector of subject numbers in the textbox under the list, such as
1, 3, 4, 5, 10, 12, 14, 15, 17<BR>
15. Analogously for level2 you may add a group/level name (Male2) with subjects (17, 18, 21, 23, 25, 29, 30).
16. Analogously for level3 you may add a group/level name (Females) with subjects (2, 6, 7, 8, 9, 11, 13, 16, 19, 20, 22, 24, 26, 27, 28).<BR>
17. Back to the "Figure 1:dFNC Stats GUI" window, you may get more results lowering the "Enter threshold in windows (Max...) from 10 to 1.<BR>
18. Click [Calculate] button.<BR>
19. Within a minute the MATLAB command window will dispaly Stats completed. Please see results file , followed by the mat file containing your statistics. This file may be loaded using the MATLAB load command.<BR>
20. To see a report with the ANOVA statistics you may click the HTML Report button from the menu at the top. You can see the MATLAB working on a report a few minutes. Finally a report will be generated in your internet browser (or you may neet to look for an html folder with a web page with ANOVA results).
