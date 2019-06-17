function kepler_elements(a, e, i, L, w, omega, mu, t)
    theta = omega*pi/180;
    g = (w - omega)*pi/180;
    i = i*pi/180;

    Tita = [ cos(theta) , -sin(theta), 0 ;
             sin(theta) ,  cos(theta), 0 ;
                 0     ,      0    , 1  ] ;
             
             
       I = [ cos(i) , 0  , -sin(i) ;
               0    , 1  ,    0    ;
             sin(i) , 0  ,  cos(i)  ] ; 
             
             
       G = [ cos(g) , -sin(g) , 0 ;
             sin(g) ,  cos(g) , 0 ;
               0    ,    0    , 1  ] ;
        
    Q = Tita*I*G;
    
    gama = 1 + mu;
    n = sqrt(gama/a^3);
    to = ((w - L)/n)*pi/180;
    l = n*(-t*2*pi - to);
    u = l + e*sin(l + e*sin(l + e*sin(l)));
    
    r = Q*a*[cos(u)-e ; sin(u)*sqrt(1-e^2) ; 0 ];
    
    v = Q*[-sin(u);cos(u)*sqrt(1-e^2);0]*a*n/(1-e*cos(u));
   
    disp('Coordinates (r)')
    disp(num2str(r))
    disp(['|r| = ', num2str(norm(r))])
   
    disp('Speed (v)')
    disp(num2str(v))
    disp(['|v| = ', num2str(norm(v))])
end