function elements(a, e, i, L, w, omega, mu, t)
     i = i*pi/180;
     n = sqrt(1/a^3);
     to = ((w - L)/n)*pi/180;    
    
     gamma = 1 + mu;
     L = mu*sqrt(gamma*a)
     G = L*sqrt(1-e^2)
     bigTheta = G*cos(i)
     l = n*(t*2*pi - to)
     g = (w - omega)*pi/180
     theta = omega*pi/180
     H =- mu*gamma/(2*a)
   
     %First Poincare System of equations
     P11 = sqrt(a)
     P12 = (L-G)/(mu*sqrt(gamma))
     P13 = (G-bigTheta)/(mu*sqrt(gamma))
     P14 = L*pi/180
     P15 = -g-theta
     P16 = -theta
     
     %Second Poincare System of equations
     P21 = P11
     P22 = sqrt(2*(L-G))*cos(g+theta)
     P23 = sqrt(2*(G-bigTheta))*cos(theta)
     P24 = P14
     P25 = -sqrt(2*(L-G))*sin(g+theta)
     P26 = -sqrt(2*(G-bigTheta))*sin(theta)
end