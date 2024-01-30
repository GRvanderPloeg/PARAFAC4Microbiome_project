function Factors=pickModel(As,Bs,Cs,choiceNumber)
if nargin<4
    error('All inputs must be given.')
end

Factors = cell(1,3);

Factors{1} = As(:,:,choiceNumber);
Factors{2} = Bs(:,:,choiceNumber);
Factors{3} = Cs(:,:,choiceNumber);