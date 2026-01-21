#State Guided ICA Cheat Sheet

Before running these State Guided measures and statistics you have to have run a blind or Neuromark ICA to perform following staps on.

1. Start the GIFT GUI.
2. Click [Dynamic Functional Connectivity...]. (see image below)
 
3. Click [Temporal dFNC (ICA)].
4. Click [Setup/Run Analysis].
5. A file dialog appears, where you select the parameter file from your first (main) ICA run.
6. A file dialog appears, where you select the output directory.
7. A window appears where you use the [ + ] button to select the components to analyze
(see the domains a and b in the attached image, which you may choose as SC, AU, SM, VI, CC,… etc… instead ).
8. Click [Run].
9. Wait until the MATLAB Command Window shows “analysis complete”.
10. Click [Post-processing].
11. Select the {prefix}_dfnc.mat file.
12. In the Post-processing dFNC parameters window (see image below):
 
•	Keep 10 components/states, or change the number of components to the number of states you want.
•	Set the Threshold level (%) to 10 or change the number of components if needed.
•	Please do not disable or change any other settings for now, as this version is not very stable for alternative configurations (just run them too as default).
•	Click [Done].
13. Wait until the Command Window shows “Done”.
14. Click [Display] 
•	Now select the dFNC param file , looking like {prefix}_dfnc.mat.
•	Now you get the Display dFNC window popping up (image below)
 
o	Select “State Guided ICA” from the list and click [Display] again (see image)
15. Report windows should now pop up (see below), 
  
•	including FNCs/States (see image), that you can click [->] to see the next report image. 
o	On the windows top you can see prior (x), where x is the component/state number
o	Over the FNC you see Pos Y seconds (Z%), Neg…, where Y shows how many seconds this FNC is active (and Z shows same number in percent) 
o	If you chose 10 components/states you will have to click [->] ten times to see the frequency and drew times (see image below)
 
o	On each window top you see Pos stats in one window and Neg stats in second window
	Stats are how frequently each of the 10 states are and the mean dwell time spent in each of the 10 states
	Analogous for the negative states
	Final report page is the Reference Guided Spatia; dFNC Transitions window (see below)
 
 
