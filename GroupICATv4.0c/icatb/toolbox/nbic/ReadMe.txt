1. Intorduction:

N-BiC: is a biclustering method for simultaneously clustering two interacting variables. A bicluster is a two-dimensional submatrix
It clusters a 2D matrix where each dimension represents a variable. Each variable has a set of examples. Variables could be subejcts, neuro components
- It requires a sorting method to build the set of instancs for all given examples of the variable 
Since we use intersection for the merging or composing the biclusters, We need to represent variables from one dimension 
as a function of variables from other dimension (as a subset of samples from other dimension) 

The given data matrix has a dimension of #samples x #features 
We perform the feature reduction to represent each feature as a (reduced) subset of samples 
Another rationale is to represent variables from one dimension as a function of variables from other dimension 
Here, we sort the components (features) as a subset of subejects based on a predefiend heuristic. 
We presented several approaches to perform the sorting in the referred paper.
However, the sorting is more intuitively carried out based on the objectives of a study - more study dependent.

     
- 'Search_BIC' is used to explore the biclusters for a given set of variables 
- The validator is integrated with the 'Search_BIC' subroutine
 

2. The subroutines:
 
The architecture is distributed into three major scripts 
   i.  NBiC: the main function. To run the N-BiC method we need to run this script with appropriate inputs settings. This script will take care of everything 
             to run the process. Finally, it returns a list of biclusters 
			 
   ii. Search_BIC: It create all possible subsets of given set of variables and creates biclusters among them. It's a recursive subroutine that uses
                   a modified DFS technique to traverse all the branches of search tree. It also uses early abondaning to kill some inadiquate branches   
                     
   iii. BiC_validation: Validate the newly formed biclusters by checking the similarity between early reported BICs. It uses F1 similarity index for contrlling the replication 
   
   iv. Feature_reduction: sort the features for reduced number of samples
   

3. Running the code base:

- We can treat the data as a 2D matrix. Each dimension of the matrix represents a variable i.e., test subjects, neuro components of a ICA loading matrix. 
  N-BiC biclusters both dimensions simultaneously 
- We need a sorting method maybe domain specific to represent/define variables from one dimension as a subset of variables from other dimension  
- Send the IDS of sorted variables and the DATA array (of these sorted variables) to the NBiC function     


Remarks:

- Further filtering might requires on the final list of biclusters based on the slection of input arguments
- Be mindful to check the frequency of the biclusters across the permutations   


Reference:
M. A. Rahaman et al., "N-BiC: A Method for Multi-Component and Symptom Biclustering of Structural MRI Data: Application to Schizophrenia," 
in IEEE Transactions on Biomedical Engineering, vol. 67, no. 1, pp. 110-121, Jan. 2020, doi: 10.1109/TBME.2019.2908815.