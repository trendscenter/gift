# N-BIC 
### Clustering Software (MATLAB)
![TReNDS](https://trendscenter.org/wp-content/uploads/2019/06/background_eeg_1.jpg)
### Table of Contents
1. [Introduction](#secIntro)
2. [Manual](#secMan)
3. [Screen Shots](#secScreen)
4. [Version History](#secVerHist)
---
### Introduction <a name="secIntro"></a>
NBIC is a clustering algorithm NBiC toolbox, based on the 2020 publication "N-BiC: A Method for Multi-Component and Symptom Biclustering of Structural MRI Data: Application to Schizophrenia" (Md Abdur Rahaman , Jessica A. Turner, Cota Navin Gupta, Srinivas Rachakonda, Jiayu Chen , Jingyu Liu , Theo G. M. van Erp, Steven Potkin, Judith Ford, Daniel Mathalon, Hyo Jong Lee, Wenhao Jiang, Bryon A. Mueller, Ole Andreassen, Ingrid Agartz, Scott R. Sponheim , Andrew R. Mayer, Julia Stephen , Rex E. Jung, Jose Canive, Juan Bustillo, and Vince D. Calhoun). This toolbox works on MATLAB versions greater than R2008a. N-BiC is a biclustering method for simultaneously clustering two interacting variables. A bicluster is a two-dimensional submatrix. It clusters a 2D matrix where each dimension represents a variable. Each variable has a set of examples. Variables could be subejcts, neuro components. It uses a sorting method based on user defiend heuristic and data requirements to select a subset of instances for defining a variable. The work horse of N-BiC algorthm is graph-search followed by intersection operation for the merging/composing the biclusters. Independent component analysis and blind source separation of group (and single subject) functional magnetic resonance imaging data. GIFT works on MATLAB R2008a and higher. Many ICA algorithms were generously contributedby Dr. Andrzej Cichocki. These are also available in Dr. Cichocki's ICALAB toolbox. For any question or comments please contact Md Abdur Rahaman (mrahaman1@gsu.edu), Vince Calhoun (vcalhoun@gsu.edu) or Cyrus Eierud (ceierud@gsu.edu). 

### Manual <a name="secMan"></a>
The given data matrix has a dimension of [sample x feature] 

#### Feature reduction/sorting
We perform the feature reduction to represent each feature as a (reduced) subset of samples. Another rationale is to represent variables from one dimension as a function of variables from other dimension. Here, we sort the components (features) as a subset of subejects based on a predefiend heuristic. We presented several approaches to perform the sorting in the referred paper. However, the sorting is more intuitively carried out based on the objectives of a study - more study dependent. We need to choose which samples we want to allow for defining a feature. Let's assume we have 100 images (sampels) of traffic signals and one of the feature we extracted is color. Now, for running N-BiC, we need to sort the feature - create a subset of samples defining that feature. For instance,  to sort the feature 'color', we might select samples carrying color = 'red'.    

#### The subroutines:
 
The architecture is distributed into three major scripts 
   i.  NBiC: the main function. To run the N-BiC method we need to run this script with appropriate inputs settings. This script will take care of everything 
             to run the process. Finally, it returns a list of biclusters 
			 
   ii. Search_BIC: It create all possible subsets of given set of variables and creates biclusters among them. It's a recursive subroutine that uses
                   a modified DFS technique to traverse all the branches of search tree. It also uses early abondaning to kill some inadiquate branches   
                     
   iii. BiC_validation: Validate the newly formed biclusters by checking the similarity between early reported BICs. It uses F1 similarity index for contrlling the replication 
   
   iv. Feature_reduction: sort the features for reduced number of samples
   
#### Running the code base:

Input: We need a 2D data matrix. Each dimension of the matrix represents a variable i.e., test subjects, neuro components of a ICA loading values. N-BiC biclusters both dimensions simultaneously. 

1. Run 'Feature_reduction'- sorting method for defining one dimension as a subset of variables from other dimension  
2. Run N-BiC with sorted variables IDs and the DATA array (of these sorted variables) for biclustering the variables from both dimension.

For more details, see the example script 'running_the_code_base.m'


#### Remarks:

- Further filtering might requires on the final list of biclusters based on the slection of input arguments
- Be mindful to check the frequency of the biclusters across the permutations   
- If you want to run N-BIC in the Group ICA interface it may be done under source-based morphometry (SBM) by:
	1. In MATLAB command window enter sbm
	2. Click [Toolbox: NBIC]
	3. Enter a parameter file (in dialog)
	4. Enter a destination directory (in dialog)
	5. Follow instructions in NBIC screen seen in Fig. 1, including adding both SBM loadings and import neuropsychological scores from csv file.
	6. Report will display in accordance with Fig. 2.

### Screen Shots <a name="secScreen"></a>

| ![GIFT](https://trendscenter.org/trends/software/gift/images/nbic1.png) |
|:--:|
| Figure 1. Setting up NBIC|

| ![GIFT](https://trendscenter.org/trends/software/gift/images/nbic2.png) |
|:--:|
| Figure 2. NBIC output figure|

### Version History<a name="secVerHist"></a>
Click the following link for the GIFT version history: [GIFT version history](https://trendscenter.org/trends/software/gift/version_history.html) 

