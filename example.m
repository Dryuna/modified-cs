
%   Example script for running MCS

%   Written by Sean Walton (512465@swansea.ac.uk) 2011 for Swansea
%   university

%   Please cite the following paper when using this code...
%   S.Walton, O.Hassan, K.Morgan and M.R.Brown "Modified cuckoo search: A
%   new gradient free optimisation algorithm" Chaos, Solitons & Fractals Vol
%   44 Issue 9, Sept 2011 pp. 710-718 DOI:10.1016/j.chaos.2011.06.004


%   I'd appreciate it if you contacted me (512465@swansea.ac.uk) if you apply the code to a
%   problem succesfully, I'm always interested in hearing about new applications 

clear all
clc

%The structure S contains all the parameters for the MCS

S.pa = 0.75;        %Fraction of eggs discarded each generation
S.maxstep = 10;     %Maximum number of steps to take in a levy flight
S.plot = 1;         %If you want the results plotted set this to 1
S.fname = 'obj';    %The function name, if this function gives a complex value the optimser considers it out of bounds
S.constrain = 1;    %Set to 1 if you want the search constrained within vardef, zero otherwise
S.A = 0.01; %   Maximum distance a cuckoo can travel in one step as fraction of search space diagonal

%   Tips for choosing the above paramters :-
%       1) S.A multiplied by S.maxstep is roughly the furthest an egg will move in a
%       generation, so is dependent on the size of your search space.  You
%       don't really want this to be larger than the diagonal of your
%       search space
%       2) When writing the function S.fname make sure it accepts a vector
%       size x(1,NoDimensions), and returns a complex value when x is
%       either invalid or outside the search space
%       3) You will get better performance when each dimension of your
%       search space is simular sized, try and formulate your function with
%       some kind of normalisation so this is the case
%       
%       


%The matrix vardef defines the upper and lower bounds of the initial set of
%nests, the MCS uses this to set boundaries on the plots and LHC uses it to
%generate initial eggs

NoDim = 10;

vardef(1,1:NoDim) = 100;
vardef(2,1:NoDim) = -100;

NoNests = 100;

NestI = LHC(vardef,NoNests); %Generates initial set of eggs

NoGen = 1000;

%Run optimiser
[p,F,pg] = MCS(NoGen, NestI, vardef, S);

%The optimum position is then pg

