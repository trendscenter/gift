%*********************************************************************
                                                                                
%********wchange.m: tracks size and direction of weight changes ******
                                                                                
%*********************************************************************
                                                                                
                                                                                
                                                                                
function [change,delta,angle]=wchange(w,oldw,olddelta)
                                                                                
  [M,N]=size(w); delta=reshape(oldw-w,1,M*N);
                                                                                
  change=delta*delta';
                                                                                
  angle=acos((delta*olddelta')/sqrt((delta*delta')*(olddelta*olddelta')));
