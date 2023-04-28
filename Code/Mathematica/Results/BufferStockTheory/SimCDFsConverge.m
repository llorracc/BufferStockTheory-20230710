(* ::Package:: *)

(*
This program generates the figure illustrating the numerical convergence
of the CDF of c.

*)
Print["Running SimCDFsConverge.m"];

If[Not[ModelIsSolved],SolveInfHorizToToleranceAtTarget[ToleranceAtTarget=0.0001]];

NumOfPeopleToSimulate  = 100000;
NumOfPeriodsToSimulate = 80;
Simulate[NumOfPeopleToSimulate,NumOfPeriodsToSimulate,bStartVal=0.];

