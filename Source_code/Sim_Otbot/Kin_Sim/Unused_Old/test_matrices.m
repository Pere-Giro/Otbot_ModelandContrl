syms a b c d e f g h k

Rmat = [a,b,c; d,e,f; g,h,k];

Amat = [1,2,3;
        1,2,3;
        1,2,3];
    
exp1 = Rmat*Amat;

disp(exp1)

Amat2 = sym(zeros(3,3));

for i=1:max(length(Amat))
    Amat2(:,i) = Rmat*Amat(:,i); % Rotation
end

disp(Amat2)