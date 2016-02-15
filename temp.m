clear all;
clc;

G1 = 7;% octal 7 corresponds to binary 111 n1 = m1 + m0 + m-1 
G2 = 3;% octal 3 corresponds to binary 011 n1 = m0  + m-1 
G3 = 5;% octal 5 corresponds to binary 101 n1 = m1  + m-1 
constLen = 3;   % Constraint length 

% Create the trellis that represents the convolutional code
convCodeTrellis = poly2trellis(constLen, [ G1 G2 G3 ]);
uncodedWord = [1];
codedWord1 = convenc(uncodedWord, convCodeTrellis) 
codedWord1 = [1 0 1];

uncodedWord = [1 0 0 0];
codedWord2 = convenc(uncodedWord, convCodeTrellis) 