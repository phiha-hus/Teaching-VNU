function test_interp(index)
 % Test interpolation scheme 
  
  switch index    
     case 1
      x = [1 -4 0]
      y = [3 13 -23]      
      case 2
        x = [1 1.5 0 2]
        y = [3 13/4 3 5/3]   
        case 3
        x = [1 2 3 4]
        y = [2 1 6 47]   
        case 4        
        x = [-2 -1 0 1 2 3]
        y = [1 4 11 16 13 -4]  
  end
    
    v = div_diff(x,y)
    w = 3;
    tic
    value = newton_interp(x,y,w)
    toc 
%
%    tic    
%    value = neville_interp(x,y,w)
%    toc
    