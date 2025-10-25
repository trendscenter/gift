# ANOVA Implemented Into MANCOVAN 10/24/2025

A 1-way (3-level) ANOVA has now been implemented in the GIFT GUI. You may find the following group-analysis example helpful: (https://github.com/trendscenter/gift/tree/master/GroupICAT/demo/gift_group_analysis).<BR>
<BR>
In order to run this you need to install the latest version (v4.0.6.13) of GIFT found at:<BR>
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
CONTINUING WITH ANOVA MODEL<BR>
3. In GIFT, click the [Toolboxes/Mancovan] drop down box<BR>
4. Click [Create Design Matrix] and then a the "Select ICA/Mancovan Parameter..." window appears<BR>
5. In the "Select..." window, browse to the output directory where you saved your parameter file from step 1, (if you downloaded the demodata from our group-analysis page (above) you pick neuromark_ica_parameter_info.mat from the demo_input3neuromark/gift_out directory.<BR>
Demo data link is: https://trends-public-website-fileshare.s3.amazonaws.com/public_website_files/software/gift/data/demo_input3neuromark.zip<BR>
6. Choose an output directory for ANOVA results.<BR>
7. On the Setup Design page:<BR>
* Select Disign Criteria: 1-way, x-level ANOVA<BR>
* When prompted for "ANOVA Levels," enter 3 (or 2 if testing equivalence to our two-sample t-test in the demo/gift_group_analysis link, above).<BR>
8. Additional windows appears where you select the subject IDs for each group (levels) and in our 2 level anova we typed "males" for group 1 and then pasted in the subject ids for the males (1, 3, 4, 5, 10, 12, 14, 15, 17, 18, 21, 23, 25, 29, 30) in the textbox above the [OK] button. After the paste  we just pressed [OK]<BR>
9. Then we typed "females" (2nd level) for group 2 and then pasted in 2, 6, 7, 8, 9, 11, 13, 16, 19, 20, 22, 24, 26, 27, 28 and then pressed [OK]<BR>
10. Click [Create]<BR>
11. (Back to the mancovan Toolbox) we then press [Setup Features]<BR>
12. Pick the mancovan file that should appear in the file list (using our data it was neuromark_mancovan.mat)<BR>
13. Select Features, which may be: Spatial Maps, Timecoruses Spectra or FNC correlations (pick one)<BR>
<BR>14. For "Add components" you click the [+] button to select the components that make your domains (for each domain, which most often is used for constrained ICA), but even if you do not use domains you at least need to pick one domain with all the components you want to report.
<BR>15. A Select component window appears and you may type "my components of interest" and then select the components in the list on the right side. To select multiple components you may use the control-key (ctrl) on your keyboard, while selecting the components of interest in the list.
<BR>16. At "Enter P-value Significance Threshold, enter a p-value you consider is a good threshold
<BR>17. Enter the TR of your analysis
<BR>18. Click "Autoselect No. Of Components For Each Feature Using MDL" to select it and then tick it again to deselect it to get the "Enter No. Of Components For Each Feature in a vector" and enter "2" in the textbox above the [Run]-button
<BR>19. Click [Run]
<BR>20. (Back to the mancovan Toolbox), press [Run Mancovan]
<BR>21. Window "About Remove nuisance" ops up and in this case click [No] (yes is mostly for multi site studies)
<BR>22. Again pick your mancovan file that should appear in the file list (using our data it was neuromark_mancovan.mat)
<BR>23. Click [Display]
<BR>24. Again pick your mancovan file that should appear in the file list (using our data it was neuromark_mancovan.mat)
<BR>25. for the label "Select results to display", pick "Univariate Results" from the dropdown
<BR>26. Select or leave settings for "Leave T-Threshold (Tmap)", "Select Image values...", "Low and high frequency...", "Display connectogram..."
<BR>27. For "Threshold Criteria.." you may want to choose "none" in case you have low statistical power
<BR>28. For "P-Threshold" you can set to 0.05 depending on your model. If you do not get any figures out you may need to set this p-value to 1 as the figures will not generate if nothing is significant and then at least you can see the insifnificant results in case you have a p=0.06 result (as a sanity check).
<BR>29. Click [Display]

You should now see ANOVA results for your selected components (from step 14) and FNC correlations
