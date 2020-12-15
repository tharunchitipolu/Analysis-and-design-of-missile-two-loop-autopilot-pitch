
        
numerator = [15.278  42.49];
denominator = [1 20.72 49.9];
sys = tf(numerator,denominator);
pzplot(sys);
rlocus(sys)

